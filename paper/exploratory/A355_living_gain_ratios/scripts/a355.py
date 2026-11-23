"""A355 (registered in `paper/docs/registered_experiments.md`): are the GTEx organ gain ratios present in living tissue?

Living cohorts with three or more tissues per person; per pair (a / b) Delta = log_amp_rel(a) - log_amp_rel(b) from
`clocksync.Dataset.fit`, both targets read against one common host (the other compartments), against the GTEx log ratio
of the same organs (`presentation/output/organs/couplings.csv`).
  GSE119717  omental / subcutaneous fat, host: the two perivascular depots      GTEx +0.465
  GSE39491   gastric cardia / squamous oesophagus, host: Barrett mucosa         GTEx -0.293
             (one host compartment: cross-moment ratio, `cross_ratio`, as `Dataset.fit` needs two anchors)
  Wu 2020    exposed (forearm, cheek) / protected (buttock) skin, host: the other layer   GTEx -0.211
A clock gene absent from a cohort is dropped from its waveform (A255 rule; GSE119717: NR1D1 and NPAS2 filtered out of
the deposited matrix by the submitters; subjects 50 and 51: duplicated subcutaneous arrays dropped). Arrays as pseudo-counts 1e6 exp(ln 2 (y - compartment median)), depth 1e8 (as `compute_chen_windows.py`).

  --stage power   gates, swap nulls (1000 per cohort; no Delta computed), SD, MDE, decidability -> output/{gates,power}.csv
  --stage score   (power note dated) Delta, swap p, bootstrap 90 %, pooled T, H355.1-2, D355.1-2 -> output/score_*.csv

    python a355.py --stage power|score [--jobs 4]
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(EXP / "A254_paired_lesions/scripts"))
sys.path.insert(0, str(EXP / "A251_skin_sites/scripts"))
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, N_SWAP, N_BOOT, N_SETS, MDE_K, N_GATE = 355, 1000, 500, 200, 2.49, 1000
WF = cs.gtex()
CORE = WF.drop(["PER1", "PER2", "NR1D1"])
GTEX_PAIRS = {"GSE119717": ("Adipose - Visceral (Omentum)", "Adipose - Subcutaneous"),
              "GSE39491": ("Stomach", "Esophagus - Mucosa"),
              "wu2020": ("Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)")}
COHORTS = list(GTEX_PAIRS)
warnings.filterwarnings("ignore")


# ---- data -------------------------------------------------------------------------------------------------------
def gse119717():
    """samples x symbols (log2) and sheet (person, compartment)."""
    import prep_lesions_tissues as prep
    s = pd.read_csv(EXP / "scouts_2026-11-12/output/paired_lesion_GSE119717_samples.csv", dtype=str)
    comp = s.status.map({"lesion (dilated AAA wall PVAT)": "pvat_dilated", "non-lesion (non-dilated neck PVAT)":
                         "pvat_neck", "subcutaneous abdominal fat": "subcutaneous", "omental visceral fat": "omental"})
    e, _ = prep.array_expr("GSE119717", "paired_lesion", "GSE119717_series_matrix.txt.gz", "GPL10558.annot.gz",
                           list(s.matrix_column))
    sheet = pd.DataFrame({"sample": s.matrix_column.to_numpy(), "person": s.patient_id.to_numpy(),
                          "compartment": comp.to_numpy()})
    # subjects 50 and 51 carry two arrays labelled subcutaneous each (scout note): ambiguous labels, both dropped
    sheet = sheet[~sheet.duplicated(["person", "compartment"], keep=False)].reset_index(drop=True)
    return e.loc[sheet["sample"]], sheet


def gse39491():
    prep = EXP / "A254_paired_lesions/output/prep"
    e = pd.read_parquet(prep / "GSE39491_expr.parquet")
    s = pd.read_csv(prep / "GSE39491_sheet.csv", dtype=str)
    sheet = pd.DataFrame({"sample": s["sample"], "person": s.person,
                          "compartment": s.condition.map({"lesion": "barrett", "nonlesion_cardia": "cardia",
                                                          "nonlesion_squamous": "squamous"})})
    return e.loc[sheet["sample"]], sheet


def wu2020():
    import compute_a251 as a251
    out, rows = [], []
    for layer, (m, s) in a251.load().items():
        w = s[s.arm != "timed_men"]
        out.append(m.loc[w.index])
        rows.append(pd.DataFrame({"sample": w.index, "person": w.subject_id.astype(str),
                                  "compartment": layer + "_" + w.site.astype(str)}))
    return pd.concat(out), pd.concat(rows, ignore_index=True)


def waveform(e, core=False):
    """The GTEx waveform on the clock genes present in the cohort (A255 rule: an absent gene is dropped there);
    core: without PER1, PER2, NR1D1."""
    wf = WF.subset([g for g in WF.genes if g in e.columns])
    return wf.drop([g for g in ("PER1", "PER2", "NR1D1") if g in wf.genes]) if core else wf


def pseudobulk(e, sheet, wf=WF):
    """Long table for `Dataset`: pseudo-counts 1e6 exp(ln2 (y - compartment median)), depth 1e8."""
    g = [x for x in wf.genes if x in e.columns]
    y = e.loc[sheet["sample"], g].astype(float).set_axis(sheet.index)
    y = y - y.groupby(sheet.compartment).transform("median")
    pb = (1e6 * np.exp(np.log(2) * y)).assign(donor=sheet.person.to_numpy(), compartment=sheet.compartment.to_numpy(),
                                              depth=1e8)
    return pb.dropna()


LOADERS = {"GSE119717": gse119717, "GSE39491": gse39491, "wu2020": wu2020}


# ---- estimate ---------------------------------------------------------------------------------------------------
def design(cohort):
    """list of (host compartments, targets a, targets b): Delta = mean log_amp(a) - mean log_amp(b) over the fits."""
    if cohort == "GSE119717":
        return [(["pvat_dilated", "pvat_neck"], ["omental"], ["subcutaneous"])]
    if cohort == "GSE39491":
        return [(["barrett"], ["cardia"], ["squamous"])]
    return [([f"{h}_{s}" for s in ("forearm", "buttock", "cheek")], [f"{t}_forearm", f"{t}_cheek"], [f"{t}_buttock"])
            for t, h in (("epidermis", "dermis"), ("dermis", "epidermis"))]


def cross_ratio(pb, host, a, b, wf):
    """One host compartment (`Dataset.fit` needs >= 2 anchors): gain ratio as log |sum z_a conj z_h| - log |sum z_b
    conj z_h| over persons with all three, states on the log pseudo-counts centred within compartment."""
    y = np.log(pb[list(wf.genes)].astype(float)).assign(donor=pb.donor.to_numpy(), compartment=pb.compartment.to_numpy())
    L = P.projector(y[list(wf.genes)], wf, counts=False)
    z = {}
    for k in (host, a, b):
        t = y[y.compartment == k].set_index("donor")[list(wf.genes)]
        t = t - t.mean()
        bb = t.to_numpy() @ L.T
        z[k] = pd.Series(bb[:, 0] + 1j * bb[:, 1], index=t.index)
    ids = z[host].index.intersection(z[a].index).intersection(z[b].index)
    h = np.conj(z[host][ids].to_numpy())
    return float(np.log(abs(np.sum(z[a][ids].to_numpy() * h))) - np.log(abs(np.sum(z[b][ids].to_numpy() * h))))


def delta(pb, cohort, wf=WF):
    d = []
    for host, a, b in design(cohort):
        if len(host) == 1:
            d.append(cross_ratio(pb, host[0], a[0], b[0], wf))
            continue
        sub = pb[pb.compartment.isin(host + a + b)]
        tab = cs.Dataset(sub, anchors=host, waveform=wf, min_donors=10).fit(n_boot=0, gate=False, targets=a + b).table
        d.append(tab.loc[a, "log_amp_rel"].mean() - tab.loc[b, "log_amp_rel"].mean())
    return float(np.mean(d))


def swapped(pb, cohort, rng):
    """Per person, the target labels exchanged (two targets: probability 1/2; Wu: the three sites permuted per layer)."""
    pb = pb.copy()
    groups = [a + b for _, a, b in design(cohort)]
    for persons_comps in groups:
        m = pb.compartment.isin(persons_comps)
        for _, i in pb[m].groupby("donor").groups.items():
            comps = pb.loc[i, "compartment"].to_numpy()
            pb.loc[i, "compartment"] = comps[rng.permutation(len(comps))]
    return pb


def swap_null(pb, cohort, k, n, wf=WF):
    return np.array([delta(swapped(pb, cohort, np.random.default_rng([SEED, k, j])), cohort, wf) for j in range(n)])


def boot(pb, cohort, k, n, wf=WF):
    ids = pb.donor.unique()
    out = []
    for j in range(n):
        rng = np.random.default_rng([SEED, 50 + k, j])
        pick = rng.choice(ids, len(ids))
        b = pd.concat([pb[pb.donor == d].assign(donor=f"{d}#{r}") for r, d in enumerate(pick)])
        try:
            out.append(delta(b, cohort, wf))
        except Exception:  # noqa: BLE001  a replicate without enough donors in a compartment
            out.append(np.nan)
    return np.array(out)


def gtex_ratio():
    d = pd.read_csv(PAPER / "presentation/output/organs/couplings.csv").set_index("tissue").log_amp_rel
    return {c: float(d[a] - d[b]) for c, (a, b) in GTEX_PAIRS.items()}


# ---- stages -----------------------------------------------------------------------------------------------------
def power(jobs):
    G, data, gates = gtex_ratio(), {}, []
    for k, c in enumerate(COHORTS):
        e, sheet = LOADERS[c]()
        wf = waveform(e)
        pb = pseudobulk(e, sheet, wf)
        data[c] = (pb, wf)
        g = cs.Dataset(pb, anchors=None, waveform=wf, min_donors=10).gate(n_perm=N_GATE, seed=SEED + k)
        n = pb.groupby("compartment").donor.nunique().to_dict()
        gates.append(dict(cohort=c, genes=len(wf.genes), p_pairing=g.p_pairing, clock_share=g.clock_share, share_null=g.share_null_mean,
                          compartments=str(n), passed=g.p_pairing < 0.05))
        print(gates[-1], flush=True)
    gates = pd.DataFrame(gates)
    gates.to_csv(OUT / "gates.csv", index=False)
    ok = list(gates[gates.passed].cohort)
    nulls = dict(zip(ok, Parallel(n_jobs=jobs)(delayed(swap_null)(*data[c][:1], c, COHORTS.index(c), N_SWAP, data[c][1]) for c in ok)))
    np.savez(OUT / "swap_nulls.npz", **nulls)
    rows = []
    for c in ok:
        sd = float(np.nanstd(nulls[c]))
        rows.append(dict(test=c, gtex_log_ratio=G[c], null_sd=sd, mde=MDE_K * sd, decidable=MDE_K * sd <= abs(G[c]),
                         n_valid=int(np.isfinite(nulls[c]).sum())))
    w = {c: 1 / np.nanstd(nulls[c]) ** 2 for c in ok}
    T0 = sum(w[c] * np.sign(G[c]) * nulls[c] for c in ok) / sum(w.values())
    target = sum(w[c] * abs(G[c]) for c in ok) / sum(w.values())
    rows.append(dict(test="pooled T", gtex_log_ratio=target, null_sd=float(np.nanstd(T0)),
                     mde=MDE_K * float(np.nanstd(T0)), decidable=MDE_K * float(np.nanstd(T0)) <= target,
                     n_valid=int(np.isfinite(T0).sum())))
    pw = pd.DataFrame(rows)
    pw.to_csv(OUT / "power.csv", index=False)
    print(pw.round(3).to_string(index=False))


def var_ratio(e, sheet, a, b, wf):
    """Total clock-plane log variance ratio a / b (lists of compartments, averaged) on log expression (ln), persons
    with both, states centred within compartment; one projector over the pair."""
    g = [x for x in wf.genes if x in e.columns]
    y = e.loc[sheet["sample"], g].astype(float).set_axis(sheet.index) * np.log(2)
    L = P.projector(y, wf if len(g) == len(wf.genes) else wf.subset(g), counts=False)

    def tr(k):
        t = y[sheet.compartment == k].to_numpy()
        bb = (t - t.mean(0)) @ L.T
        return float(np.sum(bb ** 2) / (len(t) - 1))
    return float(np.mean([np.log(tr(x)) for x in a]) - np.mean([np.log(tr(x)) for x in b]))


def companions(e, sheet, cohort, k):
    """D355.2: clock and 200 expression-matched sets, total variance ratio per target pair (averaged over the design)."""
    lm = (e.astype(float) * np.log(2)).mean()
    lm = lm[~lm.index.duplicated()].dropna()
    wf = waveform(e)
    sets = P.pseudo_waveforms(lm, N_SETS, np.random.default_rng([SEED, 90, k]), wf, "matched")
    out = []
    for w in [wf] + list(sets):
        out.append(np.mean([var_ratio(e, sheet, a, b, w) for _, a, b in design(cohort)]))
    return out[0], np.array(out[1:])


def score(jobs):
    G = gtex_ratio()
    pw = pd.read_csv(OUT / "power.csv").set_index("test")
    nulls = dict(np.load(OUT / "swap_nulls.npz"))
    ok = [c for c in COHORTS if c in nulls]
    rows, desc, data = [], [], {}
    for c in ok:
        e, sheet = LOADERS[c]()
        wf, wfc = waveform(e), waveform(e, core=True)
        pb = pseudobulk(e, sheet, wf)
        data[c] = (e, sheet, wf, wfc, pb)
    obs = {c: delta(data[c][4], c, data[c][2]) for c in ok}
    bs = dict(zip(ok, Parallel(n_jobs=jobs)(delayed(boot)(data[c][4], c, COHORTS.index(c), N_BOOT, data[c][2])
                                            for c in ok)))
    np.savez(OUT / "boot.npz", **bs)
    for c in ok:
        s = np.sign(G[c])
        p = float((1 + np.sum(s * nulls[c] >= s * obs[c])) / (1 + len(nulls[c])))
        lo, hi = np.nanquantile(bs[c], [0.05, 0.95])
        rows.append(dict(test=c, delta=obs[c], gtex_log_ratio=G[c], fold_living=float(np.exp(obs[c])),
                         fold_gtex=float(np.exp(G[c])), p_one_sided=p, lo90=lo, hi90=hi,
                         gtex_inside=bool(lo <= G[c] <= hi), decidable=bool(pw.loc[c, "decidable"]),
                         genes=len(data[c][2].genes), boot_valid=int(np.isfinite(bs[c]).sum())))
    t = pd.DataFrame(rows)
    dec = t[t.decidable]
    m = len(dec)
    t["p_holm"] = np.nan
    if m:
        order = dec.p_one_sided.sort_values().index
        adj = np.maximum.accumulate([min(1, (m - i) * t.loc[j, "p_one_sided"]) for i, j in enumerate(order)])
        t.loc[order, "p_holm"] = adj
    t["H355_2"] = np.where(t.decidable, np.where((t.p_holm < 0.05) & t.gtex_inside, "PASS", "FAIL"), "not scored")
    w = {c: 1 / float(pw.loc[c, "null_sd"]) ** 2 for c in ok}
    T = sum(w[c] * np.sign(G[c]) * obs[c] for c in ok) / sum(w.values())
    T0 = sum(w[c] * np.sign(G[c]) * nulls[c] for c in ok) / sum(w.values())
    pT = float((1 + np.sum(T0 >= T)) / (1 + len(T0)))
    Tb = sum(w[c] * np.sign(G[c]) * bs[c] for c in ok) / sum(w.values())
    pooled = dict(test="pooled T (H355.1)", delta=float(T), gtex_log_ratio=float(pw.loc["pooled T", "gtex_log_ratio"]),
                  p_one_sided=pT, lo90=float(np.nanquantile(Tb, 0.05)), hi90=float(np.nanquantile(Tb, 0.95)),
                  decidable=bool(pw.loc["pooled T", "decidable"]),
                  H355_1=("PASS" if pT < 0.05 else "FAIL") if pw.loc["pooled T", "decidable"] else "not scored",
                  weights=str({c: round(w[c] / sum(w.values()), 3) for c in ok}))
    t.to_csv(OUT / "score_pairs.csv", index=False)
    pd.DataFrame([pooled]).to_csv(OUT / "score_pooled.csv", index=False)
    for c in ok:
        e, sheet, wf, wfc, pb = data[c]
        k = COHORTS.index(c)
        pbc = pseudobulk(e, sheet, wfc)
        dc = delta(pbc, c, wfc)
        bc = boot(pbc, c, k, 200, wfc)
        vc, vs = companions(e, sheet, c, k)
        desc.append(dict(cohort=c, core_genes=len(wfc.genes), core_delta=dc, core_lo90=float(np.nanquantile(bc, 0.05)),
                         core_hi90=float(np.nanquantile(bc, 0.95)), clock_log_var_ratio=vc,
                         matched_median=float(np.median(vs)), matched_q05=float(np.quantile(vs, 0.05)),
                         matched_q95=float(np.quantile(vs, 0.95)), clock_pct_among_matched=float(np.mean(vs < vc)),
                         missing_genes=str([g for g in WF.genes if g not in wf.genes])))
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.round(3).to_string(index=False))
    print(pd.DataFrame([pooled]).round(3).to_string(index=False))
    print(pd.DataFrame(desc).round(3).to_string(index=False))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args(argv)
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "power":
        power(a.jobs)
    else:
        score(a.jobs)


if __name__ == "__main__":
    main()
