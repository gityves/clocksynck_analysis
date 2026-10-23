"""Species-specific vs human clock reference (A41): is the human GTEx waveform close to optimal in another species?

Pre-registered as A41 (2026-10-16; docs/PLAN_internal_desync.md (archive); the A36 / A40 / A34 results were seen). Baboon waveform:
`waveform.fit_waveform` on the baboon timed data (log((count + 0.5) / depth) of the 14 clock genes, ZT, group = tissue,
63 tissues).
  Readout 1 (descriptive): baboon vs human waveform gene by gene: phase after the common shift, R of the differences,
     largest |difference|, amplitude Spearman (expected if conserved: all genes within +-2 h, amplitude Spearman > 0.5).
  Readout 2 (primary): the A40 windows of 5 animals (and 6) refitted with a baboon waveform fitted only on the animals
     outside the window; host untimed; compared paired with the human-reference values. Human reference "close to
     optimal" if the median paired difference in tissue-offset Fisher-Lee r (baboon - human) < 0.10, flips correct with
     the baboon reference <= human + 1, and the median |Spearman| differs by < 0.1. Truth: full-day timed offsets with the
     human waveform; sensitivity: with the full-day baboon waveform.
  Readout 3 (secondary, in sample): the full day refitted with the baboon waveform (ordering, tissue offsets).
  Readout 4 (secondary): sheep (A34 data, untimed) refitted with the full-day baboon waveform: gate, tissue couplings
     (Spearman with the human-reference couplings), offsets agreement between references. The A34 primary (sheep offsets
     against human gains d_t) is not ported (needs the legacy sheep-human tissue matching and human gain table).
Windows, estimator and truth are those of windows/compute_baboon_windows.py (A40).

    python compute_baboon_reference.py [--subsample N] [--seed S] [--jobs 4] [--n-perm N]
`--subsample N` runs the first N start positions of each window size.
Outputs (paper/results/extended_data/output/gate/): baboon_reference_{genes,windows,summary,readouts,waveform_full}.csv
Legacy: legacy/desync/validation/baboon/scripts/{a41_reference,a41_sheep}.py -> output/a41_*.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats as sps  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import fit_waveform  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _gate_shared as sh  # noqa: E402

FOLDER = "gate"
W = sh.W
G = list(cs.gtex().genes)
SIZES = (5, 6)


def baboon_waveform(pb, donors, animals, name):
    s = pb[pb.donor.isin(animals)]
    y = np.log((s[G].to_numpy(float) + 0.5) / s.depth.to_numpy(float)[:, None])
    return fit_waveform(y, donors.zt.reindex(s.donor).to_numpy(float) % 24, G, group=s.compartment.to_numpy(), name=name)


def compare_wf(wf, ref):
    d = np.asarray(wf.psi) - np.asarray(ref.psi)
    m = np.exp(1j * d).mean()
    res = np.angle(np.exp(1j * (d - np.angle(m)))) / W
    g = pd.DataFrame(dict(gene=list(wf.genes), amp=wf.amp, psi_h=(np.asarray(wf.psi) / W) % 24, ref_amp=ref.amp,
                          ref_psi_h=(np.asarray(ref.psi) / W) % 24, resid_h=res))
    i = int(np.argmax(np.abs(res)))
    return dict(shift_h=(np.angle(m) / W) % 24, R=abs(m), max_abs_resid_h=abs(res[i]), max_gene=wf.genes[i],
                n_within_2h=int((np.abs(res) <= 2).sum()), amp_spearman=sps.spearmanr(wf.amp, ref.amp).statistic), g


def agree(a_rad, b_rad, n_perm, seed=0):
    r = cs.timed_agreement(a_rad, b_rad, n_perm=n_perm, seed=seed)
    return dict(r_fl=r.r_fl, p_fl=r.p_fl, r_js=r.r_js, p_js=r.p_js, R=r.R, offset_h=r.offset_h,
                median_abs_err_h=r.median_abs_err_h, n=int(r.n))


def offsets_vs(off_h: pd.Series, truth: pd.Series, n_perm):
    ok = off_h.notna() & truth.notna()
    return agree(off_h[ok].to_numpy(float) * W, truth[ok].to_numpy(float), n_perm) if ok.sum() >= 6 else dict(n=int(ok.sum()))


def fit_window(pb, kept, timed, tb, animals, ztu, label, arm, wf, n_perm):
    limit_threads(1)
    row = dict(window=label, arm=arm, m=len(animals), start_zt=ztu[0])
    sub = pb[pb.donor.isin(animals)].reset_index(drop=True)
    try:
        ds = cs.Dataset(sub, anchors=kept, waveform=wf, min_donors=5, name=label)
        e, ft = sh.T.host_estimate(ds, kept)
    except Exception as ex:
        return row | dict(status=f"infeasible: {type(ex).__name__}: {str(ex)[:80]}")
    t = e["tangent"]
    rs = sps.spearmanr(t.to_numpy(), [dict(zip(animals, ztu))[a] for a in t.index]).statistic
    row |= dict(status="ok", kappa=e["kappa"], r_axis=e["r"], spearman=rs, abs_spearman=abs(rs), flip_ok=bool(rs > 0))
    off = ft.table.offset_h.reindex(kept)
    row |= {f"loo_{k}": v for k, v in offsets_vs(off, timed, n_perm).items()}
    row |= {f"sens_{k}": v for k, v in offsets_vs(off, tb, n_perm).items()}
    return row


def sheep_readout(run, n_perm, wf_baboon, seed):
    spb, _ = sh.s02.built("sheep")
    wf_h = cs.gtex()
    cols = [c for c in ("donor", "compartment", "depth", "stratum") if c in spb] + G
    out, off, cpl = {}, {}, {}
    for lab, wf in (("human_ref", wf_h), ("baboon_ref", wf_baboon)):
        ds = cs.Dataset(spb[cols], anchors=sorted(spb.compartment.unique()), waveform=wf, mixed_strata=True, min_donors=5, name="sheep")
        g = ds.gate(n_perm=n_perm, seed=seed)
        out[lab] = dict(clock_share=g.clock_share, share_null=g.share_null_mean, p_pairing=g.p_pairing,
                        p_acrophase=g.p_acrophase, passes=g.passes)
        t = ds.fit(n_boot=0, gate=False).table
        off[lab], cpl[lab] = t.offset_h, t.coupling
    rows = [dict(readout=f"4a gate {k}", human_ref=out["human_ref"][k], baboon_ref=out["baboon_ref"][k]) for k in out["human_ref"]]
    ok = off["human_ref"].notna() & off["baboon_ref"].notna()
    a, b = off["baboon_ref"][ok].to_numpy() * W, off["human_ref"][ok].to_numpy() * W
    e = np.exp(1j * (a - b))
    res = np.angle(e * np.conj(e.mean())) / W
    rows += [dict(readout="4b coupling Spearman baboon-ref vs human-ref", baboon_ref=sps.spearmanr(cpl["baboon_ref"][ok], cpl["human_ref"][ok]).statistic, n=int(ok.sum())),
             dict(readout="4b coupling median", human_ref=cpl["human_ref"].median(), baboon_ref=cpl["baboon_ref"].median()),
             dict(readout="4b offsets baboon-ref vs human-ref: JS r", baboon_ref=cs.stats.js_corr(a, b), n=int(ok.sum())),
             dict(readout="4b offsets baboon-ref vs human-ref: Fisher-Lee r", baboon_ref=cs.stats.fisher_lee(a, b)),
             dict(readout="4b offsets baboon-ref vs human-ref: median |diff| after shift (h)", baboon_ref=float(np.median(np.abs(res))))]
    return rows


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=None, help="permutations (default 2000; 100 with --subsample)")
    ap.add_argument("--skip-sheep", action="store_true")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = args.n_perm if args.n_perm is not None else (100 if run["quick"] else 2000)
    pb, donors = sh.s02.baboon()
    kept = sorted(pb.compartment.unique())
    animals = sorted(pb.donor.unique())
    hum = cs.gtex()
    wf_full = baboon_waveform(pb, donors, animals, "baboon_full")
    r1, genes = compare_wf(wf_full, hum)
    genes.to_csv(out / "baboon_reference_genes.csv", index=False)
    wf_full.table().to_csv(out / "baboon_reference_waveform_full.csv")
    rows = [dict(readout="1 baboon full-day waveform vs human", **r1)]
    print(genes.round(3).to_string(index=False))
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r1.items()})
    timed = sh.s02.timed_tissue_gains(pb, donors).timed_phase.reindex(kept)
    tb = sh.s02.timed_tissue_gains(pb, donors, wf=wf_full).timed_phase.reindex(kept)
    tt = agree(tb.to_numpy(), timed.to_numpy(), n_perm)
    rows.append(dict(readout="truth check: timed phase baboon-wf vs human-wf (63 tissues)", **tt))
    # readout 2
    jobs, wfs = [], {}
    n_starts = 12 if args.subsample is None else min(12, max(2, args.subsample // 6))
    for m in SIZES:
        for s in range(n_starts):
            an = [animals[(s + i) % 12] for i in range(m)]
            out_an = [a for a in animals if a not in an]
            ztu = [2 * s + 2 * i for i in range(m)]
            lab = f"m{m}_ZT{2 * s:02d}-ZT{(2 * s + 2 * (m - 1)) % 24:02d}"
            wfs[lab] = baboon_waveform(pb, donors, out_an, f"baboon_out_{lab}")
            jobs += [(an, ztu, lab, "human", hum), (an, ztu, lab, "baboon", wfs[lab])]
    wrows = []
    for lab, wf in wfs.items():
        c, _ = compare_wf(wf, hum)
        c2, _ = compare_wf(wf, wf_full)
        wrows.append(dict(window=lab, **c, R_vs_full=c2["R"], max_abs_resid_vs_full_h=c2["max_abs_resid_h"]))
    WW = pd.DataFrame(wrows)
    res = Parallel(n_jobs=args.jobs)(delayed(fit_window)(pb, kept, timed, tb, *j, n_perm) for j in jobs)
    D = pd.DataFrame(res).merge(WW[["window", "R", "max_abs_resid_h", "amp_spearman"]].rename(
        columns=dict(R="wf_R", max_abs_resid_h="wf_max_resid_h", amp_spearman="wf_amp_spearman")), on="window")
    D.to_csv(out / "baboon_reference_windows.csv", index=False)
    S = []
    for m in SIZES:
        h = D[(D.arm == "human") & (D.m == m)].set_index("window")
        b = D[(D.arm == "baboon") & (D.m == m)].set_index("window")
        ok = [w for w in b.index if b.status.get(w) == "ok" and h.status.get(w) == "ok"]
        if not ok:
            continue
        dlt = b.loo_r_fl[ok] - h.loo_r_fl[ok]
        S.append(dict(m=m, n_paired=len(ok), loo_r_fl_human_med=h.loo_r_fl[ok].median(), loo_r_fl_baboon_med=b.loo_r_fl[ok].median(),
                      delta_r_fl_med=dlt.median(), delta_r_fl_min=dlt.min(), delta_r_fl_max=dlt.max(), n_baboon_better=int((dlt > 0).sum()),
                      loo_mae_human_med=h.loo_median_abs_err_h[ok].median(), loo_mae_baboon_med=b.loo_median_abs_err_h[ok].median(),
                      flips_human=int(h.flip_ok[ok].astype(bool).sum()), flips_baboon=int(b.flip_ok[ok].astype(bool).sum()),
                      abs_spearman_human_med=h.abs_spearman[ok].median(), abs_spearman_baboon_med=b.abs_spearman[ok].median(),
                      sens_r_fl_human_med=h.sens_r_fl[ok].median(), sens_r_fl_baboon_med=b.sens_r_fl[ok].median(),
                      sens_delta_r_fl_med=(b.sens_r_fl[ok] - h.sens_r_fl[ok]).median()))
    S = pd.DataFrame(S)
    S.to_csv(out / "baboon_reference_summary.csv", index=False)
    pd.set_option("display.width", 220)
    print(S.round(3).T.to_string())
    if len(S) and (S.m == 5).any():
        p = S[S.m == 5].iloc[0]
        c1, c2 = p.delta_r_fl_med < 0.10, p.flips_baboon <= p.flips_human + 1
        c3 = abs(p.abs_spearman_baboon_med - p.abs_spearman_human_med) < 0.1
        print(f"PRIMARY (m=5): median paired delta FL r {p.delta_r_fl_med:.3f} (< 0.10: {c1}); flips {p.flips_baboon} vs "
              f"{p.flips_human} (<= +1: {c2}); |Spearman| {p.abs_spearman_baboon_med:.3f} vs {p.abs_spearman_human_med:.3f} "
              f"(< 0.1: {c3}) -> {'human reference CLOSE TO OPTIMAL' if c1 and c2 and c3 else 'human reference COSTS accuracy'}")
    # readout 3: full day, baboon waveform (in sample)
    ds = cs.Dataset(pb, anchors=kept, waveform=wf_full, name="baboon_full_bwf")
    ft = ds.fit(n_boot=0, gate=False)
    zt = donors.zt.astype(float)
    r = ds.timed(zt, n_perm=n_perm)
    rows.append(dict(readout="3 full day baboon-wf: host phase vs ZT", n=int(r.n), R=r.R, p_R=r.p_R, r_fl=r.r_fl, p_fl=r.p_fl,
                     r_js=r.r_js, p_js=r.p_js, offset_h=r.offset_h, median_abs_err_h=r.median_abs_err_h))
    off_b = ft.table.offset_h.reindex(kept)
    rows.append(dict(readout="3 full day baboon-wf: offsets vs timed (human-wf truth)", **offsets_vs(off_b, timed, n_perm)))
    rows.append(dict(readout="3 full day baboon-wf: offsets vs timed (baboon-wf truth)", **offsets_vs(off_b, tb, n_perm)))
    dh = cs.Dataset(pb, anchors=kept, name="baboon_full_human")
    th = dh.fit(n_boot=0, gate=False).table
    rows.append(dict(readout="3 A36 human-wf offsets vs timed (human-wf truth)", **offsets_vs(th.offset_h.reindex(kept), timed, n_perm)))
    rows.append(dict(readout="3 A36 human-wf offsets vs timed (baboon-wf truth)", **offsets_vs(th.offset_h.reindex(kept), tb, n_perm)))
    rows.append(dict(readout="3 full day: offsets baboon-wf vs human-wf", **agree(off_b.to_numpy() * W, th.offset_h.reindex(kept).to_numpy() * W, n_perm)))
    rows.append(dict(readout="3 full day: coupling Spearman baboon-wf vs human-wf",
                     r_fl=sps.spearmanr(ft.table.coupling, th.coupling).statistic, n=len(kept)))
    if not args.skip_sheep:
        rows += sheep_readout(run, run["n_perm"], wf_full, args.seed)
    R = pd.DataFrame(rows)
    R.to_csv(out / "baboon_reference_readouts.csv", index=False)
    print(R.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
