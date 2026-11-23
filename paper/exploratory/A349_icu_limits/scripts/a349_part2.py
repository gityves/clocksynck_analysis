"""A349 part 2 (registered, descriptive): host spread, host-referenced log z2 fast vs ventilated post-mortem, relative gain.

GTEx fast deaths (Hardy 1-2) vs the 79 ventilated post-mortem donors with a living brain (A293 part 1 group), 12 body
tissues of A109 as targets, host anchored on the 15 A164 tissues (as A306), `clocksync.Dataset`; 14 genes and core loop.
  D349.1  host spread V per group, separate fits (`Dataset.fit`; V of the targets, median over targets, sigma_f trace),
          and V per group from the pooled-host contrast.
  D349.2  per tissue log z2 ratio fast / ventilated and d* at the common spread V0 (`Dataset.contrast`), 90 % intervals
          from a 500-draw donor bootstrap within groups (host and V0 fixed, as the contrast's own bootstrap, percentile
          intervals); pooled and matched on sex, age (caliper 5 y) and ischaemic time (A335 rule, `match_ratio`).
  D349.3  relative organ gain = per tissue log z2 ratio minus its median over tissues; Spearman of the per-tissue loss
          with the tissue's fast-death z2.
Aggregates only (GTEx controlled access): no donor-level output.

    python a349_part2.py [--check]
Outputs (`../output/`): part2_spread.csv, part2_tissues.csv, part2_summary.csv
"""
import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import a349_common as C

sys.path.insert(0, str(C.EXP / "A231_age_coupling" / "scripts"))
from a349_common import B, H, S, gtex
import clocksync as cs
from clocksync._core import host as A
from clocksync import stats as ST
from clocksync.config import limit_threads
import matched_label as ML
import _shared as SH

TRIO = ["PER1", "PER2", "NR1D1"]
LEVELS = ("ventilated", "fast")  # case = fast: ratio fast / ventilated


def datasets(ids):
    pb, _ = gtex.load(tissues=list(SH.TISSUES), hardy=(0, 1, 2), genes="clock", depth="library_size", require_tod=True)
    pb = pb[pb.donor.isin(ids)].reset_index(drop=True)
    anc = [a for a in SH.TISSUES if a in set(pb.compartment)]
    return {"14 genes": cs.Dataset(pb, anchors=anc, name="a349_14"),
            "core": cs.Dataset(pb, anchors=anc, name="a349_core", drop_genes=TRIO)}


def boot(ds, lab, targets, n, seed):
    """Percentile bootstrap (donors within groups, host and V0 fixed) of log z2 ratio and d*: {tissue: (lz2 draws, d* draws)}."""
    lab = lab.astype(str)
    d = ds.data
    d = d.take(np.flatnonzero(np.isin(d.resection, lab.index)))
    need = list(targets)
    Hs = A._hosts(d, [a for a in ds.anchors if a in set(d.compartment)], need)
    idx = {k: pd.Index(np.intersect1d(d.donors(k), Hs[k].resections)) for k in need}
    V0 = {k: A.on_host(d, k, Hs[k], idx[k].to_numpy())["V"] for k in need}
    rng = np.random.default_rng(seed)
    out = {}
    ratio = lambda a, b: np.log(a / b) if (a or 0) > 0 and (b or 0) > 0 else np.nan  # noqa: E731
    for k in need:
        lk = lab.reindex(idx[k]).dropna()
        don = {g: lk.index[lk == g].to_numpy() for g in LEVELS}
        rows_k = {r: np.flatnonzero((d.compartment == k) & (d.resection == r)) for r in np.concatenate(list(don.values()))}
        reps = []
        for _ in range(n):
            r = []
            for g in LEVELS:
                e = A.on_host_drawn(d, k, Hs[k], rng.choice(don[g], len(don[g])), rows_k)
                s = e["z2"] * V0[k]
                r.append((e["z2"], e["psi"], s / (s + 2 * e["psi"]) if np.isfinite(s) else np.nan))
            reps.append([ratio(r[1][0], r[0][0]), r[1][2] - r[0][2]])
        out[k] = np.array(reps, float)
    return out


def run(ds, lab, gs, design, n_boot, seed):
    tab = ds.contrast(lab, levels=LEVELS, targets=H.BODY, n_perm=0, n_boot=0, blocks="none").table
    b = boot(ds, lab, [t for t in H.BODY if t in tab.index and tab.loc[t, "status"] == "ok"], n_boot, seed)
    rows = []
    for t in H.BODY:
        if t not in b:
            continue
        x = b[t]
        r = tab.loc[t]
        rows.append(dict(gene_set=gs, design=design, tissue=t, n_fast=int(r["donors_fast"]), n_vent=int(r["donors_ventilated"]),
                         log_z2_ratio=r.log_z2_ratio, lz2_lo90=np.nanquantile(x[:, 0], 0.05), lz2_hi90=np.nanquantile(x[:, 0], 0.95),
                         d_star=r.d_star, dstar_lo90=np.nanquantile(x[:, 1], 0.05), dstar_hi90=np.nanquantile(x[:, 1], 0.95),
                         V0=r.V0, V_fast=r.V_fast, V_ventilated=r.V_ventilated, z2_fast=r.z2_fast, z2_ventilated=r.z2_ventilated,
                         log_psi_ratio=r.log_psi_ratio, c_fast=r.c_fast, c_ventilated=r.c_ventilated,
                         n_boot_ok=int(np.isfinite(x[:, 0]).sum())))
    t = pd.DataFrame(rows)
    t["relative_gain"] = t.log_z2_ratio - t.log_z2_ratio.median()
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    out = C.OUT / ("check" if a.check else ".")
    out.mkdir(parents=True, exist_ok=True)
    nb = 3 if a.check else 500
    T, fast, pm = C.tissues()
    ids_f, ids_p = T.donors[fast].astype(str), T.donors[pm].astype(str)
    dss = datasets(set(ids_f) | set(ids_p))
    print(f"data ({time.time() - t0:.0f}s): fast {len(ids_f)}, ventilated {len(ids_p)}", flush=True)
    # matching: A335 rule (sex exact, age caliper 5 y, nearest age then ischaemic time, up to 3 controls per case)
    ph = gtex.phenotypes(B.PH_COLS + ["AGE", "SEX"]).reindex(T.donors)
    attrs = S.sample_attrs(T)
    med = lambda col: np.array([np.nanmedian([attrs[t][col].to_numpy(float)[T.R[t][d]] for t in H.BODY if T.R[t][d] >= 0]  # noqa: E731
                                             or [np.nan]) for d in range(len(T.donors))])
    q = pd.DataFrame(dict(age=pd.to_numeric(ph.AGE, errors="coerce").to_numpy(),
                          male=(pd.to_numeric(ph.SEX, errors="coerce") == 1).astype(int).to_numpy(), hardy=0,
                          ischaemic_h=med("ischemic_time") / 60), index=np.arange(len(T.donors)))
    sets = ML.match_ratio(q, pd.Series(np.isin(q.index, pm), index=q.index), pd.Series(np.isin(q.index, fast), index=q.index))
    cd = np.array([s[0] for s in sets])
    ctl = np.array(sorted({x for s in sets for x in s[1:]}))
    print(f"matched: {len(cd)} ventilated, {len(ctl)} fast controls", flush=True)
    labs = {"pooled": pd.Series(["fast"] * len(ids_f) + ["ventilated"] * len(ids_p), index=list(ids_f) + list(ids_p)),
            "matched": pd.Series(["fast"] * len(ctl) + ["ventilated"] * len(cd),
                                 index=list(T.donors[ctl].astype(str)) + list(T.donors[cd].astype(str)))}
    spread, tabs, summ = [], [], []
    for gs, ds in dss.items():
        for g, ids in (("fast", ids_f), ("ventilated", ids_p)):
            f = ds.subset(donors=ids).fit(n_boot=0, gate=False)
            v = f.table.reindex(H.BODY).V
            spread.append(dict(gene_set=gs, group=g, n_donors=len(ids), V_median_targets=float(v.median()), V_min=float(v.min()),
                               V_max=float(v.max()), sigma_f_trace=float(np.trace(f.sigma_f)), **{f"V_{t}": v[t] for t in H.BODY}))
        for design, lab in labs.items():
            t = run(ds, lab, gs, design, nb, [349, 2, C.SEEDS[gs], design == "pooled"])
            tabs.append(t)
            rho = spearmanr(t.log_z2_ratio, t.z2_fast)
            summ.append(dict(gene_set=gs, design=design, n_tissues=len(t), median_log_z2_ratio=float(t.log_z2_ratio.median()),
                             n_tissues_ratio_gt0=int((t.log_z2_ratio > 0).sum()), median_V_fast=float(t.V_fast.median()),
                             median_V_ventilated=float(t.V_ventilated.median()), median_V0=float(t.V0.median()),
                             median_log_V_ratio=float(np.log(t.V_fast / t.V_ventilated).median()),
                             median_d_star=float(t.d_star.median()), median_c_fast=float(t.c_fast.median()),
                             median_c_ventilated=float(t.c_ventilated.median()),
                             relative_gain_min=float(t.relative_gain.min()), relative_gain_max=float(t.relative_gain.max()),
                             spearman_loss_vs_z2_fast=float(rho[0]), spearman_p=float(rho[1])))
            print(f"{gs} / {design} done ({time.time() - t0:.0f}s)", flush=True)
    pd.DataFrame(spread).to_csv(out / "part2_spread.csv", index=False)
    pd.concat(tabs).to_csv(out / "part2_tissues.csv", index=False)
    pd.DataFrame(summ).to_csv(out / "part2_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(spread).iloc[:, :8].round(4).to_string(index=False))
    print(pd.DataFrame(summ).round(3).T.to_string())
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
