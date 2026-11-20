"""A336 part 1, Figure 6: A332 D1 (psoriasis loudness: lesion over non-lesion clock-plane variance ratio r), descriptive.

Original: `A332_descriptive/scripts/d1_skin_loudness.py` (imported, not edited). Rows: GSE154200 psoriasis (A332_D1) and the four-series
A115 pooled psoriasis row (A332_D1_A115pooled). Same data, seed 332 (`SEED + tag` for the cohort fits, [332, tag, 4|5] companions,
10,000 swaps, 200 sets); only the waveform changes. Companions: the original 14-gene draws restricted to the core positions
(A329 rule). The entry is descriptive (no verdict): core_status 'descriptive'.

    python p1_a332.py [--check]     # --check: tiny run, assertions only, nothing printed about values nor written
Output: ../output/part1/a332.csv.
"""
import argparse
import importlib.util
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE.parent))
spec = importlib.util.spec_from_file_location("p1_d1", PAPER / "exploratory/A332_descriptive/scripts/d1_skin_loudness.py")
D1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(D1)
import a336_common as H  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

ENTRY = "a332"
SCRIPT = "paper/exploratory/A336_two_readouts/scripts/p1_a332.py"
SER = pd.read_csv(PAPER / "exploratory/A332_descriptive/output/d1_series.csv").set_index("group")
SEED = D1.SEED
UNIT = dict(unit="unit", condition="condition", levels=("nonlesional", "lesional"), cohort="cohort")


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def fit(samples, counts, wf_read, tag):
    """P.paired cohort table and the gated keys for waveform wf_read (genes absent dropped), as the original Group."""
    present = [g for g in wf_read.genes if g in samples.columns]
    res = P.paired(samples, waveform=wf_read.subset(present), counts=counts, n_perm=1, n_boot=1000, seed=SEED + tag, **UNIT)
    keys = list(res.cohorts[res.cohorts.gated & np.isfinite(res.cohorts.r) & (res.cohorts.se > 0)].cohort)
    return res, keys


def pooled_r(res, keys):
    t = res.cohorts.set_index("cohort").loc[keys]
    iv = 1 / t.se.to_numpy() ** 2
    return float(np.sum(iv / iv.sum() * t.r)), iv


class Core(D1.Group):
    """The original Group with the waveform read swapped; companions drawn on the 14-gene waveform, restricted to `wf_read`'s positions."""

    def __init__(self, name, samples, counts, means, n_perm, n_sets, tag, wf_read):
        self.name, self.n_perm, self.n_sets = name, n_perm, n_sets
        present = [g for g in D1.WF.genes if g in samples.columns]
        wf14 = D1.WF.subset(present)
        self.wf = wf_read.subset([g for g in wf_read.genes if g in samples.columns])
        res, self.keys = fit(samples, counts, wf_read, tag)
        self.cohorts, self.res, self.data, self.comp = res.cohorts, res, res.data, {}
        if not self.keys:
            return
        for k, x in samples.groupby("cohort"):
            if k in self.keys:
                case, ref = P._split(x, "unit", "condition", ("nonlesional", "lesional"), counts)
                c = counts if isinstance(counts, bool) else bool((x[counts] == "counts").all())
                for kind, t in (("matched", 4), ("responsive", 5)):
                    rng = np.random.default_rng([SEED, tag, t])
                    ws = P.pseudo_waveforms(means[k], n_sets, rng, wf14, kind)
                    keep = H.core(wf14)
                    self.comp.setdefault(kind, {})[k] = np.array([P.log_var_ratio(*P.cohort_states(case, ref, c, H.restrict(w, wf14, keep))) for w in ws])


def build(n_perm, n_sets, mode, wf_read=None, quick=False):
    """{label: Group} for GSE154200 psoriasis and the four A115 series (+ pooled); mode '14' = the original class, 'core' = Core."""
    A = D1.a183()
    s = A.load()
    means = A.gene_means(s)
    mk = (lambda *a, **k: D1.Group(*a)) if mode == "14" else (lambda *a, **k: Core(*a, wf_read))
    out = {}
    g = mk("GSE154200 psoriasis", s[s.diagnosis == "psoriasis"], True, means, n_perm, n_sets, 10)
    out["GSE154200 psoriasis"] = g
    counts_of, mm, series = {}, {}, {}
    for i, k in enumerate(D1.A115):
        x = D1.a115_samples(k, counts_of)
        nl = x[x.condition == "nonlesional"]
        mm[k] = D1.gene_means(nl.assign(depth=nl.depth.to_numpy(float)) if counts_of[k] else nl, counts_of[k])
        series[k] = mk(f"A115 {k}", x, counts_of[k], mm, n_perm, n_sets, 20 + i)
    ok = [k for k in D1.A115 if series[k].keys]
    pool = D1.Group.__new__(D1.Group)
    pool.n_perm, pool.n_sets = n_perm, n_sets
    pool.cohorts = pd.concat([series[k].cohorts for k in ok], ignore_index=True)
    pool.data = {k: series[k].data[k] for k in ok}
    pool.comp = {kind: {k: series[k].comp[kind][k] for k in ok} for kind in ("matched", "responsive")}
    pool.keys = ok
    out["A115 pooled (inverse variance)"] = pool
    return out


def rows_of(G):
    return {lab: (g.row(lab) if g.keys else None) for lab, g in G.items() if lab in ("GSE154200 psoriasis", "A115 pooled (inverse variance)")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(2)
    t0 = time.time()
    n_perm, n_sets, nb = (30, 3, 10) if chk else (10_000, 200, 200)
    G14 = build(n_perm, n_sets, "14")
    R14 = rows_of(G14)
    log("stage A done")
    bad = {}
    for lab, r in R14.items():
        for k in ("r", "se", "p_swap_two_sided", "matched_q05", "matched_q95", "responsive_q05", "responsive_q95", "matched_percentile", "responsive_percentile"):
            if not np.isclose(r[k], SER.loc[lab, k], rtol=0, atol=1e-6):
                bad[f"{lab}.{k}"] = (r[k], SER.loc[lab, k])
    if not chk:
        log("reproduced" if not bad else f"NOT reproduced {bad}")
    if not chk and bad:
        H.write_rows(ENTRY, [dict(id=i, figure_panel="Fig 6 (psoriasis loudness)", statistic=f"r, {lab}", direction="+1", value_14_stored=SER.loc[lab, "r"], value_14=R14[lab]["r"],
                                  reproduced=False, notes=f"14 genes did not reproduce d1_series.csv; no core value computed. {bad}", script=SCRIPT)
                             for i, lab in (("A332_D1", "GSE154200 psoriasis"), ("A332_D1_A115pooled", "A115 pooled (inverse variance)"))])
        return
    WF = H.WF
    core, trio = H.core(WF), H.trio(WF)
    Gc = build(n_perm, n_sets, "core", core)
    Rc = rows_of(Gc)
    # unit bootstrap (pairs within cohort), shared by both waveforms
    rng = np.random.default_rng([336, 1])
    d = {lab: [] for lab in G14}
    pooled_lab = ("GSE154200 psoriasis", "A115 pooled (inverse variance)")
    iv = {lab: {} for lab in pooled_lab}
    for lab in pooled_lab:
        for name, G in (("14", G14), ("core", Gc)):
            t = G[lab].cohorts.set_index("cohort").loc[G[lab].keys]
            iv[lab][name] = (list(t.index), 1 / t.se.to_numpy() ** 2)
    for _ in range(nb):
        picks = {}
        for lab in pooled_lab:
            for name, G in (("14", G14), ("core", Gc)):
                for k in G[lab].keys:
                    if (lab, k) not in picks:
                        n = len(G[lab].data[k][0])
                        picks[(lab, k)] = rng.integers(0, n, n)
        for lab in pooled_lab:
            vals = {}
            for name, G in (("14", G14), ("core", Gc)):
                ks, w = iv[lab][name]
                rb = [P.log_var_ratio(*[x[picks[(lab, k)]] for x in G[lab].data[k]]) for k in ks]
                vals[name] = float(np.nansum(w / w.sum() * rb))
            d[lab].append(vals["14"] - vals["core"])
    # gene-drop null (r pooled by inverse variance over the gated cohorts; no companions, no nulls)
    A = D1.a183()
    s = A.load()
    psor = s[s.diagnosis == "psoriasis"]
    a115 = {}
    cnt = {}
    for k in D1.A115:
        a115[k] = D1.a115_samples(k, cnt)
    dl, tl = H.drops(WF)
    wfs = dict(dl)
    v1, v2 = {}, {}
    for lab, w in (dl[:3] if chk else dl):
        try:
            res, ks = fit(psor, True, w, 10)
            v1[lab] = pooled_r(res, ks)[0] if ks else np.nan
            rs, ivs = {}, {}
            tot_r, tot_iv = 0.0, 0.0
            for i, k in enumerate(D1.A115):
                resk, kk = fit(a115[k], cnt[k], w, 20 + i)
                if kk:
                    t = resk.cohorts.set_index("cohort").loc[kk]
                    ivk = 1 / t.se.to_numpy() ** 2
                    rs[k], ivs[k] = float(t.r.iloc[0]), float(ivk[0])
            v2[lab] = float(sum(ivs[k] * rs[k] for k in rs) / sum(ivs.values())) if rs else np.nan
        except Exception:
            v1[lab] = v2[lab] = np.nan
    # trio gate (reference gate on the trio, the A183 / A332 gate)
    def tgate(samples, counts, tag):
        try:
            res, _ = fit(samples, counts, trio, tag)
            return res.cohorts.gate_p.to_numpy()
        except Exception:
            return np.array([np.nan])
    tg1 = tgate(psor, True, 10)
    if chk:
        assert tl in wfs and H.drops(WF, n=3)[1] == tl
        print("check ok")
        return
    tg2 = np.concatenate([tgate(a115[k], cnt[k], 20 + i) for i, k in enumerate(D1.A115)])
    assert np.isclose(v1[tl], Rc["GSE154200 psoriasis"]["r"]), (v1[tl], Rc["GSE154200 psoriasis"]["r"])
    assert np.isclose(v2[tl], Rc["A115 pooled (inverse variance)"]["r"]), (v2[tl], Rc["A115 pooled (inverse variance)"]["r"])
    pc = {"GSE154200 psoriasis": H.drop_percentiles(v1, tl, wfs, +1)}
    pc["A115 pooled (inverse variance)"] = H.drop_percentiles(v2, tl, wfs, +1)
    rows = []
    for i, lab, tg in (("A332_D1", "GSE154200 psoriasis", tg1), ("A332_D1_A115pooled", "A115 pooled (inverse variance)", tg2)):
        r, rc = R14[lab], Rc[lab]
        dd = np.array(d[lab])
        gate_ok = bool(np.isfinite(tg).all() and (tg < 0.05).all()) if np.isfinite(tg).any() else "not computable"
        note = (f"descriptive (no verdict); core: r {rc['r']:.4f} (amplitude ratio {rc['amplitude_ratio']:.3f}), swap p {rc['p_swap_two_sided']:.4g}, matched 90% range [{rc['matched_q05']:.3f}, {rc['matched_q95']:.3f}] "
                f"(clock at percentile {rc['matched_percentile']:.3f}), responsive [{rc['responsive_q05']:.3f}, {rc['responsive_q95']:.3f}] (percentile {rc['responsive_percentile']:.3f}), beyond companions {rc['beyond_companions']} "
                f"(14 genes: {r['beyond_companions']}); cohorts gated {rc['cohorts']}; 'beyond companions' is the entry's only criterion and carries no verdict; "
                f"trio gate: reference gate p per cohort {np.round(tg, 3).tolist()} (3 genes allow 6 acrophase arrangements: p cannot fall below about 1/6)"
                + ("; gene-drop null on the GSE154200 pooled r" if lab == "GSE154200 psoriasis" else "; gene-drop null on the inverse-variance r of the gated A115 series (fixed-effect, one cohort per series)"))
        rows.append(dict(id=i, figure_panel="Fig 6 (psoriasis loudness)", statistic=f"r (log clock-plane variance ratio lesion / non-lesion), {lab}", direction="+1",
                         value_14_stored=SER.loc[lab, "r"], value_14=r["r"], reproduced=True, value_core=rc["r"],
                         core_criterion="descriptive, no verdict ('beyond companions' = r outside the central 90% of both companion families)", core_status="descriptive",
                         diff_14_minus_core=r["r"] - rc["r"], diff_lo90=float(np.quantile(dd, .05)), diff_hi90=float(np.quantile(dd, .95)),
                         diff_draws="200 pair bootstraps within cohort (the statistic's SE bootstrap is not a draw of r), seed [336, 1], shared by both waveforms", **pc[lab],
                         trio_gate_passes=gate_ok, trio_alone_value=np.nan, reading="descriptive (no verdict)", blind="core blind", notes=note, script=SCRIPT))
    H.write_rows(ENTRY, rows)
    log("done", f"{(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
