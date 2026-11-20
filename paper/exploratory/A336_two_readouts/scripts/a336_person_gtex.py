"""A336 part 2 (registered): person loudness in GTEx deaths on the core loop.

A331's leverage (`clocksync.elasticity.radius_covariance`: residualised correlation of the radii of two disjoint
half-host states of the same donor), on the 20 random halvings of the 27 peripheral tissues of A331's power stage
(`split_hosts(..., 20, seed=[331, 1])`; fast deaths, A104 / A312 rows), median over halvings. Waveforms: 14 genes
(reproduces A331's power-stage cv2 first), core loop, 200 expression-matched and 200 responsive sets (14-gene draws,
seeds [336, 2, 1] / [336, 2, 2], restricted to the core positions for the core readout: A329 / A335 rule), every
11-gene drop (364), the trio alone (descriptive, if its gate passes). Persons-permuted null (1,000; seed [336, 2, 3]),
14-gene minus core on 1,000 shared donor bootstraps (seed [336, 2, 4]), and age / ischaemic time / RIN tertiles
(descriptive). Aggregates only (GTEx controlled access).

    python a336_person_gtex.py --stage companions [--jobs 4]   # companion leverages only (no clock value)
    python a336_person_gtex.py --stage score [--jobs 4]
    python a336_person_gtex.py --stage check                   # 2 sets, 2 drops; nothing saved or printed
Outputs (`../output/`): person_gtex_companions.csv | person_gtex.csv, person_gtex_drops.csv, person_gtex_strata.csv.
"""
import argparse
import itertools
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync import elasticity as E  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
META = ["donor", "compartment", "depth"]
WF = cs.gtex()
TRIO = ["PER1", "PER2", "NR1D1"]
CORE = WF.drop(TRIO)
CORE_IDX = [k for k, g in enumerate(WF.genes) if g in CORE.genes]
DROPS = list(itertools.combinations(WF.genes, 3))
TRIO_DROP = DROPS.index(tuple(TRIO))
N_SETS, N_PERM, N_BOOT, N_SPLITS = 200, 1000, 1000, 20
S = {}


def load():
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(a.var_names)
    a.file.close()
    pb, donors = gtex.load(min_donors=A104.MIN_DONORS, genes=genes)
    keep = ~pb.compartment.str.startswith("Brain")
    return pb[keep].reset_index(drop=True), donors


def states(i):
    """Half-host states of every halving for waveform i: list of (donors, zA, zB)."""
    warnings.filterwarnings("ignore")
    limit_threads(1)
    lab, w = S["wfs"][i]
    pb, comps = S["pb"], S["comps"]
    ds = cs.Dataset(pb[META + list(w.genes)], anchors=comps, waveform=w, name=str(lab))
    out = []
    for hA, hB in E.split_hosts(ds, comps, N_SPLITS, seed=[331, 1]):
        j = hA.index.intersection(hB.index)
        out.append((np.asarray(j), (hA.x + 1j * hA.y)[j].to_numpy(), (hB.x + 1j * hB.y)[j].to_numpy()))
    return i, out


def lev(st, keep=None, perm=None):
    """Median over halvings of the leverage (and cv2); keep: donor subset; perm: per-halving permutation of B."""
    L, C = [], []
    for h, (j, a, b) in enumerate(st):
        if keep is not None:
            m = np.isin(j, keep)
            a, b = a[m], b[m]
        if perm is not None:
            b = b[perm[h]]
        r = E.radius_covariance(a, b)
        L.append(r["leverage"])
        C.append(r["cv2"])
    return float(np.median(L)), float(np.median(C))


def boot_idx(st, rng):
    """One donor bootstrap shared by all halvings and waveforms: per halving, the row indices of the drawn donors."""
    don = np.unique(np.concatenate([j for j, _, _ in st]))
    draw = rng.choice(don, len(don))
    return [np.concatenate([np.flatnonzero(j == d) for d in draw]) for j, _, _ in st]


def run(wfs, jobs):
    S["wfs"] = wfs
    res = [None] * len(wfs)
    t0 = time.time()
    with Pool(jobs) as pool:
        for k, (i, st) in enumerate(pool.imap_unordered(states, range(len(wfs)))):
            res[i] = st
            if k % 25 == 0:
                print(f"{k + 1}/{len(wfs)} waveforms {time.time() - t0:.0f}s", flush=True)
    return res


def tertiles(v):
    v = pd.Series(v).dropna()
    r = v.rank(method="first")
    return pd.qcut(r, 3, labels=["low", "mid", "high"])


def sv(w):
    r = np.asarray(w.amp) / np.asarray(w.sd)
    return np.linalg.svd(np.c_[r * np.cos(w.psi), r * np.sin(w.psi)], compute_uv=False)


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["companions", "score", "check"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.stage == "check"
    limit_threads(1)
    t0 = time.time()
    pb, donors = load()
    comps = sorted(pb.compartment.unique())
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    lm = lm[~lm.index.duplicated()]
    del x
    n = 2 if chk else N_SETS
    sets = {kind: pseudo_waveforms(lm, n, np.random.default_rng([336, 2, tag]), WF, kind=kind)
            for kind, tag in (("matched", 1), ("responsive", 2))}
    core_of = lambda w, j: Waveform(tuple(w.genes[k] for k in CORE_IDX), CORE.amp, CORE.psi, None, None, f"core{j}")  # noqa: E731
    wfs = [((gs, kind, j), w if gs == "14 genes" else core_of(w, j))
           for kind, ws in sets.items() for j, w in enumerate(ws) for gs in ("14 genes", "core")]
    keep = sorted({g for _, w in wfs for g in w.genes} | set(WF.genes))
    S.update(pb=pb[META + keep].copy(), comps=comps)
    del pb
    print(f"{len(comps)} tissues, {S['pb'].donor.nunique()} donors, {len(keep)} genes {time.time() - t0:.0f}s", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)

    if a.stage in ("companions", "check"):
        st = run(wfs, a.jobs)
        comp = pd.DataFrame([dict(gene_set=lab[0], kind=lab[1], set=lab[2], leverage=lev(s)[0], cv2=lev(s)[1])
                             for (lab, _), s in zip(wfs, st)])
        if chk:
            print("check companions ok", flush=True)
        else:
            comp.to_csv(OUT / "person_gtex_companions.csv", index=False)
            q = comp.groupby(["gene_set", "kind"]).leverage.describe(percentiles=[0.5, 0.95])
            print(q.round(4).to_string(), flush=True)
            return

    # ---- score: 14 genes reproduced first ---------------------------------------------------------------------------
    drops = [WF.drop(list(c)) for c in DROPS]
    if chk:
        drops = drops[:1] + [drops[TRIO_DROP]]
    clock = [("14 genes", WF), ("core", CORE), ("trio", WF.subset(TRIO))]
    st = run(clock[:1], a.jobs)[0]
    l14, c14 = lev(st)
    if not chk:
        ref = float(pd.read_csv(PAPER / "exploratory/P1_person_loudness/output/power.csv", index_col=0).loc["gtex_cv2_median", "value"])
        print(f"A331 power cv2 reproduction: {c14:.6f} vs {ref:.6f}", flush=True)
        if not np.isclose(c14, ref, atol=1e-9):
            raise SystemExit("A331 power stage not reproduced: stop before the core loop")
    stc = run(clock[1:2], a.jobs)[0]
    stt = run(clock[2:], a.jobs)[0]
    std = run([(f"drop{k}", w) for k, w in enumerate(drops)], a.jobs)
    comp = pd.read_csv(OUT / "person_gtex_companions.csv") if not chk else pd.DataFrame(
        dict(gene_set=["14 genes", "core"] * 2, kind=["matched"] * 2 + ["responsive"] * 2, leverage=[0.0] * 4))
    rng = np.random.default_rng([336, 2, 3])
    perms = [[rng.permutation(len(j)) for j, _, _ in st] for _ in range(20 if chk else N_PERM)]
    rows = []
    for gs, s in (("14 genes", st), ("core", stc)):
        L, C = lev(s)
        null = np.array([lev(s, perm=p)[0] for p in perms])
        r = dict(gene_set=gs, leverage=L, cv2=C, p_perm=float((1 + np.sum(null >= L)) / (1 + len(null))),
                 null_q95=float(np.quantile(null, 0.95)), donors=len(np.unique(np.concatenate([j for j, _, _ in s]))))
        for kind in ("matched", "responsive"):
            v = comp[(comp.gene_set == gs) & (comp.kind == kind)].leverage
            r.update({f"{kind}_median": float(v.median()), f"{kind}_q95": float(v.quantile(0.95)),
                      f"pct_{kind}": float((v < L).mean() * 100)})
        r["holds"] = bool(L > r["matched_q95"] and r["p_perm"] < 0.05)
        rows.append(r)
    # paired difference on shared donor bootstraps
    rb = np.random.default_rng([336, 2, 4])
    diff = []
    for _ in range(5 if chk else N_BOOT):
        idx = boot_idx(st, rb)
        lb = [lev([(j[i], a_[i], b_[i]) for (j, a_, b_), i in zip(s, idx)])[0] for s in (st, stc)]
        diff.append(lb[0] - lb[1])
    # gene-drop null
    nd = np.array([lev(s)[0] for s in std])
    svs = np.array([sv(w) for w in drops])
    ti = len(drops) - 1 if chk else TRIO_DROP
    assert np.isclose(nd[ti], rows[1]["leverage"])
    sim = np.all(np.abs(svs / sv(CORE) - 1) <= 0.2, axis=1)
    pct = float(np.mean(nd < nd[ti]) * 100)
    # trio alone with its gate
    trio_ds = cs.Dataset(S["pb"][META + TRIO], anchors=comps, waveform=WF.subset(TRIO), name="trio")
    g = trio_ds.gate(n_perm=20 if chk else 1000, seed=336)
    gate_ok = bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05)
    rows[0].update(diff_14_minus_core=rows[0]["leverage"] - rows[1]["leverage"],
                   diff_lo90=float(np.quantile(diff, 0.05)), diff_hi90=float(np.quantile(diff, 0.95)),
                   trio_gate_p_pairing=g.p_pairing, trio_gate_p_acrophase=g.p_acrophase, trio_gate_passes=gate_ok,
                   trio_alone_leverage=lev(stt)[0] if gate_ok else np.nan)
    rows[1].update(drop_pct_all=pct, n_drops=len(nd), n_similar=int(sim.sum()),
                   drop_pct_similar=float(np.mean(nd[sim] < nd[ti]) * 100) if sim.sum() >= 20 else np.nan,
                   drop_q10=float(np.quantile(nd, 0.10)), drop_median=float(np.median(nd)), below_q10=bool(pct < 10))
    c = rows[1]
    c["reading"] = ("holds: person loudness is a clock property in deaths" if c["holds"] else
                    "falls: carried by PER1, PER2 and NR1D1" if c["below_q10"] else "falls: not decidable (precision)")
    # descriptive strata (aggregates only)
    attrs = gtex_attrs()
    don = np.unique(np.concatenate([j for j, _, _ in st]))
    cov = pd.DataFrame(index=don)
    cov["age"] = donors.reindex(don).age.to_numpy(float) if "age" in donors else np.nan
    cov["ischaemic"] = donors.reindex(don).ischaemic_h.to_numpy(float) if "ischaemic_h" in donors else np.nan
    cov["rin"] = attrs.reindex(don).to_numpy(float)
    srows = []
    for v in ("age", "ischaemic", "rin"):
        if cov[v].notna().sum() < 30:
            continue
        tt = tertiles(cov[v])
        for lvl in ("low", "mid", "high"):
            k = tt.index[tt == lvl]
            srows.append(dict(covariate=v, tertile=lvl, donors=len(k), median_value=float(cov.loc[k, v].median()),
                              leverage_14=lev(st, keep=k)[0], leverage_core=lev(stc, keep=k)[0]))
    if chk:
        print("check ok", flush=True)
        return
    pd.DataFrame(rows).to_csv(OUT / "person_gtex.csv", index=False)
    pd.DataFrame(dict(drop=["+".join(c_) for c_ in DROPS], sv1=svs[:, 0], sv2=svs[:, 1], leverage=nd)).to_csv(
        OUT / "person_gtex_drops.csv", index=False)
    pd.DataFrame(srows).to_csv(OUT / "person_gtex_strata.csv", index=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(rows).round(4).T.to_string())
    print(pd.DataFrame(srows).round(3).to_string(index=False))
    print(f"reading (registered): {c['reading']}; done {time.time() - t0:.0f}s", flush=True)


def gtex_attrs():
    """Donor median RIN over the A104 fast-death samples (first sample per donor x tissue)."""
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    o = a.obs[["subject", "tissue", "rin", "DTHHRDY"]].copy()
    a.file.close()
    o = o[o.DTHHRDY.isin(list(gtex.FAST)) & ~o.tissue.astype(str).str.startswith(("Brain", "Cells -"))]
    o["subject"] = o.subject.astype(str)
    return o.drop_duplicates(["subject", "tissue"]).groupby("subject").rin.median()


if __name__ == "__main__":
    main()
