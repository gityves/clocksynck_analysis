"""A37: living human skin, epidermis + dermis every 4 h (GSE205155; 11 subjects x 7 times): untimed relative phase and host phase against clock time.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A37; run note: (1) fails on the CI clause, (2) passes narrowly).
Hypotheses: A37 (1) the untimed dermis-minus-epidermis offset lies within 1.5 h of the timed offset and the subject-bootstrap 95 %
CI of the difference covers 0; (2) the circular correlation of the per-sample host phase with clock time is > 0.5 with p < 0.05.
Estimands: (1) untimed offset = arg of the complex cross-moment of the two layers' clock-plane states, centred within subject
(`paired.cross`; for two compartments it is the phase of the free-host loading z; |z| is not identified, as stated before
scoring); timed offset = arg of the ratio of the layers' complex loadings on clock time (local `timed_loading`, subject
intercepts). (2) `Dataset.timed` (host of the two layers, subject as stratum). Interval: subject bootstrap of (1) (closed form,
1000 draws) and of (2) (host refitted, `Dataset.resample`; 20 draws in --subsample, 100 in a full run). Null of (2): clock
times permuted over samples (not within subject, `Dataset.timed`), 2000 permutations (200 in --subsample).
Sensitivities: s1 no subject stratum, s2 1 M pseudo-reads, s3 GTEx skin waveform, s4 without the 32 h units.

    python compute_skin_truth.py [--subsample N]   (N: random subjects)
Outputs: paper/results/figure_1/output/method/skin_truth_{summary,gate}.csv, skin_truth_boot.csv
Legacy: legacy/desync/validation/skin_paired/scripts/skin_paired.py (output/summary.csv)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import cross as paired_cross  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402

FOLDER = "method"
REF, OTHER = "epidermis", "dermis"
HOURS = S.HOURS


def states(ds: cs.Dataset, donors: pd.DataFrame, centre: bool = True):
    """Per layer: clock-plane states (centred within subject), traces of their measurement covariance, unit order, times."""
    d = ds.data
    out = {}
    for k in (REF, OTHER):
        i = np.flatnonzero(d.compartment == k)
        unit = d.resection[i]
        subj = donors.subject.reindex(unit).to_numpy()
        b = d.beta[i].copy()
        for q in np.unique(subj) if centre else []:
            b[subj == q] -= b[subj == q].mean(0)
        subj_g = subj if centre else np.full(len(subj), "")
        out[k] = dict(unit=unit, subj=subj, grp=subj_g, beta=b, raw=d.beta[i], v=np.trace(d.V[i], axis1=1, axis2=2),
                      hour=donors.clock_h.reindex(unit).to_numpy(float))
    assert (out[REF]["unit"] == out[OTHER]["unit"]).all()
    return out


def offsets(st, idx=None) -> dict:
    idx = np.arange(len(st[REF]["unit"])) if idx is None else idx
    e, o = st[REF], st[OTHER]
    c = paired_cross(o["beta"][idx], o["v"][idx], e["beta"][idx], e["v"][idx])
    ze = S.timed_loading(e["raw"][idx], e["hour"][idx], e["grp"][idx])
    zd = S.timed_loading(o["raw"][idx], o["hour"][idx], o["grp"][idx])
    return dict(offset_untimed=c["offset_h"], corr_untimed=c["corr"], amp_lo=c["amp_lo"], amp_hi=c["amp_hi"],
                offset_timed=float(np.angle(zd / ze) * HOURS), amp_timed=abs(zd / ze),
                clockoff_epi=float(np.angle(ze) * HOURS % 24), clockoff_derm=float(np.angle(zd) * HOURS % 24),
                absamp_epi=abs(ze), absamp_derm=abs(zd))


def run_one(pb, donors, label, strata, waveform, n_boot, n_perm, seed, quick):
    pb = pb if strata else pb.drop(columns="stratum")
    ds = cs.Dataset(pb, anchors=[REF, OTHER], donors=donors, strata="subject" if strata else None, waveform=waveform)
    st = states(ds, donors, centre=strata)
    row = offsets(st)
    row["diff"] = float(S.cdiff(row["offset_untimed"], row["offset_timed"]))
    t = ds.timed(donors.clock_h, n_perm=n_perm, seed=seed)
    row.update(rho=t.r_js, R=t.R, offset_h=t.offset_h, mae_h=t.median_abs_err_h, p_perm=t.p_js, rho_fl=t.r_fl, p_fl=t.p_fl)
    row = dict(run=label, n_subjects=donors.subject.nunique(), n_units=len(st[REF]["unit"]), n_arrays=2 * len(st[REF]["unit"])) | row
    boot = None
    if n_boot:
        rng = np.random.default_rng(seed)
        subs = np.unique(st[REF]["subj"])
        reps = []
        for _ in range(n_boot * 10):  # closed form: 10 x more draws than the host refits
            pick = rng.choice(subs, len(subs))
            idx = np.concatenate([np.flatnonzero(st[REF]["subj"] == p) for p in pick])
            r = offsets(st, idx)
            r["diff"] = float(S.cdiff(r["offset_untimed"], r["offset_timed"]))
            for k in ("offset_untimed", "offset_timed", "diff"):
                r[k] = row[k] + float(S.cdiff(r[k], row[k]))
            reps.append(r)
        boot = pd.DataFrame(reps)
        for k in ("diff", "offset_untimed", "offset_timed"):
            row[f"{k}_lo"], row[f"{k}_hi"] = np.nanpercentile(boot[k], [2.5, 97.5])
        # host-phase statistics: subject bootstrap with the host refitted (subjects drawn with repeats are relabelled `<id>#j`)
        hb = []
        subj_units = {q: list(donors.index[donors.subject == q]) for q in np.unique(donors.subject)}
        for b in range(max(n_boot // 10, 5)):
            parts, hours = [], {}
            for j, q in enumerate(rng.choice(sorted(subj_units), len(subj_units))):
                x = pb[pb.donor.isin(subj_units[q])].copy()
                x["donor"] = x.donor + f"#{j}"
                if "stratum" in x:
                    x["stratum"] = x.stratum.astype(str) + f"#{j}"
                parts.append(x)
                hours.update({u + f"#{j}": donors.clock_h[u] for u in subj_units[q]})
            try:
                dsb = cs.Dataset(pd.concat(parts, ignore_index=True), anchors=[REF, OTHER], donors=None, waveform=waveform)
                tb = dsb.timed(pd.Series(hours), n_perm=1)
                hb.append(dict(rho=tb.r_js, R=tb.R, mae_h=tb.median_abs_err_h))
            except Exception:
                continue
        hb = pd.DataFrame(hb)
        for k in ("rho", "R", "mae_h"):
            row[f"{k}_lo"], row[f"{k}_hi"] = np.nanpercentile(hb[k], [2.5, 97.5]) if len(hb) else (np.nan, np.nan)
    return row, boot, ds


def main():
    ap = cli.parser(__doc__)
    args = ap.parse_args()
    opt = cli.setup(args)
    out = out_dir(FOLDER, opt["quick"])
    n_perm, n_boot = (200, 20) if opt["quick"] else (2000, 100)
    builder = data._builder("skin_timed_gse205155")
    pb, donors = builder.build()
    if args.subsample is not None:
        subs = np.random.default_rng(args.seed).permutation(sorted(donors.subject.unique()))[:args.subsample]
        donors = donors[donors.subject.isin(subs)]
        pb = pb[pb.donor.isin(donors.index)].reset_index(drop=True)
    pb["donor"] = pb.donor.astype(str)
    donors.index = donors.index.astype(str)
    print(f"feasible (>= 8 subjects): {donors.subject.nunique() >= 8}", flush=True)
    rows = []
    row, boot, ds = run_one(pb, donors, "primary (subject strata, 30 M, global waveform)", True, "human", n_boot, n_perm, 37, opt["quick"])
    rows.append(row)
    boot.to_csv(out / "skin_truth_boot.csv", index=False)
    g = ds.gate(n_perm=200 if opt["quick"] else 1000, seed=37)
    g.to_series().to_csv(out / "skin_truth_gate.csv", header=False)
    print(g.summary(), flush=True)
    rows.append(run_one(pb, donors, "s1 no subject stratum", False, "human", 0, n_perm, 37, opt["quick"])[0])
    pb1, d1 = builder.build(n_reads=1e6)
    pb1 = pb1[pb1.donor.isin(donors.index)].reset_index(drop=True)
    rows.append(run_one(pb1, donors, "s2 1 M pseudo-reads", True, "human", 0, n_perm, 37, opt["quick"])[0])
    rows.append(run_one(pb, donors, "s3 GTEx skin waveform", True, "Skin", 0, n_perm, 37, opt["quick"])[0])
    keep = donors.index[donors.label_h != 32]
    rows.append(run_one(pb[pb.donor.isin(keep)].reset_index(drop=True), donors.loc[keep], "s4 without 32h units", True, "human", 0, n_perm, 37, opt["quick"])[0])
    summ = pd.DataFrame(rows)
    summ.to_csv(out / "skin_truth_summary.csv", index=False)
    r = rows[0]
    p1 = abs(r["diff"]) <= 1.5 and r["diff_lo"] <= 0 <= r["diff_hi"]
    p2 = r["rho"] > 0.5 and r["p_perm"] < 0.05
    pd.set_option("display.width", 250)
    print(summ.set_index("run").T.to_string())
    print(f"(1) untimed {r['offset_untimed']:.2f} h, timed {r['offset_timed']:.2f} h, diff {r['diff']:.2f} h "
          f"[{r['diff_lo']:.2f}, {r['diff_hi']:.2f}] -> {'PASS' if p1 else 'FAIL'}")
    print(f"(2) rho {r['rho']:.3f} (perm p {r['p_perm']:.4f}), R {r['R']:.3f}, median err {r['mae_h']:.2f} h -> {'PASS' if p2 else 'FAIL'}")


if __name__ == "__main__":
    main()
