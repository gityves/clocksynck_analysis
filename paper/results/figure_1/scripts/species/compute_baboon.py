"""A36, A40, A41: baboon (GSE98965, 64 tissues around the clock) as a known-time test of the untimed host, working-hours windows,
and the species-specific clock reference. A42 (brain vs periphery in baboon and sheep) is in organs.

A36 (pre-registered): each animal is a donor (one per ZT, 12), each tissue a compartment (63 kept), the human GTEx waveform,
all tissues anchors. R1 (primary): the host phase of the 63 tissues reads the ZT of the animals (JS circular correlation > 0.5 with
permutation p < 0.05; Fisher-Lee r and R as the post hoc addendum); sensitivities with peripheral or neural / ocular anchors; leave
one animal out. R2: the untimed tissue offsets (`Dataset.fit`, leave-one-out host) against the timed tissue phase (-arg g, the
12-ZT cosinor gain on the waveform). R3: the baboon timed gain phase against the GTEx timed gain phase over fixed tissue pairs
(>= 6; r > 0, p < 0.05). Descriptive gate.
A40 (pre-registered): cyclic windows of consecutive animals (3-6, all 12 starts), readout (1): tissue offsets estimated from the
window only against the full-day timed offsets (Fisher-Lee r, median error after a common shift); primary at 8 h (5 animals):
median r > 0.5. Needs >= 5 animals per tissue (4 h and 6 h windows are infeasible with the unchanged estimator). Readouts (2) ordering
and (3) window centre use the A25 curvature estimator, not in clocksync (not ported).
A41 (pre-registered): the baboon waveform fitted on its own 12-ZT series (`fit_waveform`) against the human one; the A36 readouts
with the baboon waveform (host phase vs ZT, offsets vs the timed truth under either waveform).

    python compute_baboon.py [--subsample N] [--jobs J] [--gate]   (N: animals; checks run m = 5 windows only)
Outputs: paper/results/figure_1/output/species/a36_{animals,tissues,pairs,readouts,stability}.csv, a40_{windows,summary}.csv,
a41_{readouts,genes}.csv
Legacy: validation/baboon/scripts/{score,a40_window,a41_reference,human_gain,common}.py (A36 readout 1 stability bootstrap, A40
readouts 2-3 and the A41 window refits with curvature not ported)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats  # noqa: E402

import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "species"
W = S.W
PAIRS = {
    "ADC": "Adrenal Gland", "AOR": "Artery - Aorta", "ASC": "Colon - Transverse", "DEC": "Colon - Sigmoid",
    "HEA": "Heart - Left Ventricle", "KIC": "Kidney - Cortex", "LIV": "Liver", "LUN": "Lung",
    "MUG": "Muscle - Skeletal", "OES": "Esophagus - Muscularis", "OMF": "Adipose - Visceral (Omentum)",
    "WAS": "Adipose - Subcutaneous", "PAN": "Pancreas", "PIT": "Pituitary", "PRO": "Prostate",
    "SKI": "Skin - Not Sun Exposed (Suprapubic)", "STF": "Stomach", "TES": "Testis", "THR": "Thyroid",
    "CER": "Brain - Cerebellum", "HIP": "Brain - Hippocampus", "PUT": "Brain - Putamen (basal ganglia)",
    "AMY": "Brain - Amygdala", "PRC": "Brain - Frontal Cortex (BA9)", "SUN": "Brain - Substantia nigra",
}
BRAIN_PAIRS = {"CER", "HIP", "PUT", "AMY", "PRC", "SUN"}
NEURAL = {"AMY", "ARC", "CER", "DMH", "HAB", "HIP", "LGP", "LH", "MGP", "MMB", "OLB", "PON", "PRA", "PRC", "PUT", "PVN", "SCN", "SON",
          "SUN", "THA", "VIC", "VMH", "PIN", "PIT", "RET", "RPE", "COR", "ONH", "IRI"}  # neural / ocular tissues (anchor sensitivity)
SIZES = (3, 4, 5, 6)


def zt_rad(donors):
    return np.asarray(donors, float) * W


def agree(a, b, n_perm, seed=0):
    r = cs.timed_agreement(a, b, n_perm=n_perm, seed=seed)
    return r.to_dict()


def tissue_table(ds, pb, zt, wf, comps):
    """Untimed offsets and timed gains of every tissue (R2 inputs)."""
    t = ds.fit(n_boot=0, gate=False).table
    rows = []
    for k in comps:
        sub = pb[pb.compartment == k].set_index("donor").loc[zt.index.intersection(pb[pb.compartment == k].donor)]
        g, ng = S.timed_gain(sub, zt.reindex(sub.index).to_numpy(float), wf)
        rows.append(dict(tissue=k, neural=k in NEURAL, offset_h=t.offset_h[k], coupling=t.coupling[k], boundary=bool(t.boundary[k]),
                         amplitude=t.rel_amplitude[k], g_abs=abs(g), g_arg=np.angle(g), g_genes=ng))
    T = pd.DataFrame(rows).set_index("tissue")
    T["timed_phase"] = -T.g_arg
    T["g_phi_h"] = (T.g_arg / W) % 24
    return T


def window_job(pb, donors, comps, timed, an, label, m, start, n_perm):
    limit_threads(1)
    row = dict(window=label, m=m, start_zt=2 * start, animals=",".join(an), span_h=2 * (m - 1))
    if m < 5:
        return row | dict(status="infeasible: on_host needs >= 5 rows per tissue")
    try:
        ds = cs.Dataset(pb[pb.donor.isin(an)], anchors=comps, name=label)
        t = ds.fit(n_boot=0, gate=False).table
        z = (t.offset_h * W).reindex(comps)
        ok = z.notna()
        r = agree(z[ok].to_numpy(), timed[ok].to_numpy(), n_perm)
        return row | dict(status="ok", **{f"loo_{k}": v for k, v in r.items()})
    except Exception as ex:  # noqa: BLE001
        return row | dict(status=f"infeasible: {type(ex).__name__}: {ex}")


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--gate", action="store_true", help="descriptive gate of the 63 tissues (slow)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    out = out_dir(FOLDER, quick)
    n_perm = 1000 if quick else 10_000
    pb, donors = data.load("baboon")
    zt = donors.zt.astype(float)
    comps = sorted(pb.compartment.unique())
    wf = cs.gtex()
    ds = cs.Dataset(pb, anchors=comps, donors=donors, name="baboon")
    rows = []
    # ---- R1
    r1 = ds.timed(zt, n_perm=n_perm)
    rows.append(dict(readout="1 PRIMARY host phase vs ZT (all tissues)", **r1.to_dict()))
    for lab, anc in (("peripheral anchors", [k for k in comps if k not in NEURAL]), ("neural/ocular anchors", [k for k in comps if k in NEURAL])):
        rows.append(dict(readout=f"1 sensitivity: {lab} ({len(anc)})", **ds.timed(zt, n_perm=n_perm, anchors=anc).to_dict()))
    loao = []
    for a in (donors.index[:3] if quick else donors.index):
        sub = ds.subset(donors=[x for x in donors.index if x != a])
        r = sub.timed(zt, n_perm=200 if quick else 1000)
        loao.append(dict(left_out=a, r_js=r.r_js, r_fl=r.r_fl, R=r.R))
    pd.DataFrame(loao).to_csv(out / "a36_stability.csv", index=False)
    # ---- R2
    T = tissue_table(ds, pb, zt, wf, comps)
    T.to_csv(out / "a36_tissues.csv")
    z = T.offset_h.to_numpy() * W
    rows.append(dict(readout="2 untimed offset vs timed -arg g (all tissues)", **agree(z, T.timed_phase.to_numpy(), n_perm)))
    for lab, m in (("peripheral", ~T.neural), ("neural/ocular", T.neural)):
        rows.append(dict(readout=f"2 sensitivity: {lab}", **agree(z[m.to_numpy()], T.timed_phase[m].to_numpy(), n_perm)))
    # ---- R3: GTEx gains of the paired tissues
    hg = S.gtex_gains(sorted(set(PAIRS.values())), n_boot=0 if quick else 200, subsample=args.subsample, seed=args.seed).set_index("tissue")
    P = pd.DataFrame([dict(baboon=b, gtex=g) for b, g in PAIRS.items() if b in T.index and g in hg.index])
    P["baboon_g_phi_h"] = T.g_phi_h.reindex(P.baboon).to_numpy()
    P["baboon_g_abs"] = T.g_abs.reindex(P.baboon).to_numpy()
    P["baboon_untimed_lag_h"] = ((-T.offset_h.reindex(P.baboon) * W / W) % 24).to_numpy()
    P["human_g_phi_h"] = hg.phi_g_h.reindex(P.gtex).to_numpy()
    P["human_g_abs"] = hg.abs_g.reindex(P.gtex).to_numpy()
    P["human_n"] = hg.n.reindex(P.gtex).to_numpy()
    P.to_csv(out / "a36_pairs.csv", index=False)
    bw, hw = P.baboon_g_phi_h.to_numpy() * W, P.human_g_phi_h.to_numpy() * W
    per = ~P.baboon.isin(BRAIN_PAIRS).to_numpy()
    rows.append(dict(readout="3 cross-species: baboon timed arg g vs human arg g", **agree(bw, hw, n_perm)))
    rows.append(dict(readout="3 sensitivity: peripheral pairs", **agree(bw[per], hw[per], n_perm)))
    rows.append(dict(readout="3 sensitivity: baboon untimed (-offset) vs human arg g", **agree(P.baboon_untimed_lag_h.to_numpy() * W, hw, n_perm)))
    R = pd.DataFrame(rows)
    R.to_csv(out / "a36_readouts.csv", index=False)
    pd.set_option("display.width", 250)
    print(R.round(4).to_string(index=False))
    p1, p3 = R.iloc[0], R[R.readout.str.startswith("3 cross")].iloc[0]
    print(f"PRIMARY as written (JS r > 0.5 and p < 0.05): {'PASSES' if p1.r_js > .5 and p1.p_js < .05 else 'FAILS'} (r {p1.r_js:.3f})")
    print(f"PRIMARY with Fisher-Lee (post hoc): {'PASSES' if p1.r_fl > .5 and p1.p_fl < .05 else 'FAILS'} (r_fl {p1.r_fl:.3f})")
    print(f"Cross-species (r > 0, p < 0.05, >= 6 pairs): {'PASSES' if p3.r_js > 0 and p3.p_js < .05 and p3.n >= 6 else 'FAILS'} (r {p3.r_js:.3f})")
    if args.gate:
        print(ds.gate(n_perm=200).to_series().round(4).to_dict())

    # ---- A40 windows (readout 1)
    allan = list(donors.index)
    sizes = (5,) if quick else SIZES
    starts = range(3 if quick else 12)
    jobs = [([allan[(s + i) % 12] for i in range(m)], f"m{m}_ZT{2 * s:02d}-ZT{(2 * s + 2 * (m - 1)) % 24:02d}", m, s)
            for m in sizes for s in starts]
    timed = (T.timed_phase).reindex(comps)
    W_rows = Parallel(n_jobs=args.jobs)(delayed(window_job)(pb, donors, comps, timed, an, lab, m, s, 500 if quick else 2000)
                                        for an, lab, m, s in jobs)
    D = pd.DataFrame(W_rows)
    D.to_csv(out / "a40_windows.csv", index=False)
    summ = []
    for m, g in D.groupby("m"):
        ok = g[g.status == "ok"]
        r = dict(m=m, span_h=2 * (m - 1), n_windows=len(g), n_ok=len(ok))
        for c in ("loo_r_fl", "loo_r_js", "loo_median_abs_err_h"):
            if c in ok and ok[c].notna().any():
                r[f"{c}_med"], r[f"{c}_min"], r[f"{c}_max"] = ok[c].median(), ok[c].min(), ok[c].max()
        summ.append(r)
    summ = pd.DataFrame(summ)
    summ.to_csv(out / "a40_summary.csv", index=False)
    print(summ.round(3).to_string(index=False))
    ok5 = D[(D.m == 5) & (D.status == "ok")]
    if len(ok5):
        print(f"PRIMARY (1) m=5 median Fisher-Lee r = {ok5.loo_r_fl.median():.3f} (> 0.5): {'PASSES' if ok5.loo_r_fl.median() > .5 else 'FAILS'}")

    # ---- A41 baboon waveform
    G = list(wf.genes)
    s = pb.set_index("donor").loc[zt.index]
    s = pb[pb.donor.isin(zt.index)]
    y = np.log((s[G].to_numpy(float) + 0.5) / s.depth.to_numpy(float)[:, None])
    wfb = cs.waveform.fit_waveform(y, s.donor.map(zt).to_numpy(float) % 24, G, group=s.compartment.to_numpy(), name="baboon")
    d = np.asarray(wfb.psi) - np.asarray(wf.psi)
    m = np.exp(1j * d).mean()
    res = np.angle(np.exp(1j * (d - np.angle(m)))) / W
    i = int(np.argmax(np.abs(res)))
    genes_t = pd.DataFrame(dict(gene=G, amp=wfb.amp, psi_h=(np.asarray(wfb.psi) / W) % 24, ref_amp=wf.amp,
                                ref_psi_h=(np.asarray(wf.psi) / W) % 24, resid_h=res))
    genes_t.to_csv(out / "a41_genes.csv", index=False)
    a41 = [dict(readout="1 baboon full-day waveform vs human", shift_h=(np.angle(m) / W) % 24, R=abs(m), max_abs_resid_h=abs(res[i]),
                max_gene=G[i], n_within_2h=int((np.abs(res) <= 2).sum()), amp_spearman=stats.spearmanr(wfb.amp, wf.amp).statistic)]
    dsb = cs.Dataset(pb, anchors=comps, donors=donors, waveform=wfb, name="baboon_wf")
    r = dsb.timed(zt, n_perm=n_perm)
    a41.append(dict(readout="3 full day baboon-wf: host phase vs ZT", **r.to_dict()))
    Tb = tissue_table(dsb, pb, zt, wfb, comps)
    zb = Tb.offset_h.to_numpy() * W
    a41.append(dict(readout="3 full day baboon-wf: offsets vs timed (human-wf truth)", **agree(zb, T.timed_phase.to_numpy(), n_perm)))
    a41.append(dict(readout="3 full day baboon-wf: offsets vs timed (baboon-wf truth)", **agree(zb, Tb.timed_phase.to_numpy(), n_perm)))
    a41.append(dict(readout="3 A36 human-wf offsets vs timed (baboon-wf truth)", **agree(z, Tb.timed_phase.to_numpy(), n_perm)))
    a41.append(dict(readout="3 full day: offsets baboon-wf vs A36 human-wf", **agree(zb, z, n_perm)))
    a41.append(dict(readout="truth check: timed phase baboon-wf vs human-wf", **agree(Tb.timed_phase.to_numpy(), T.timed_phase.to_numpy(), n_perm)))
    A41 = pd.DataFrame(a41)
    A41.to_csv(out / "a41_readouts.csv", index=False)
    print(A41.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
