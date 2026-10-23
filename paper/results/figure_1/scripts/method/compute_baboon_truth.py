"""A36: baboon multi-tissue (GSE98965; 12 animals, one per ZT0-22, 63 tissues): does the untimed host read the true time?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A36; run note: the Jammalamadaka-SenGupta primary is degenerate for 12 evenly
spaced times, the Fisher-Lee correlation is the disclosed post hoc statistic). Hypotheses: A36 readout 1 (PRIMARY: circular
correlation of host phase with ZT > 0.5, p < 0.05, ZT permuted), readout 2 (untimed tissue offsets vs timed tissue phases over
63 tissues), readout 3 (cross-species: baboon timed gain phase vs human GTEx gain phase over the matched tissues, >= 6).
Estimand: R, circular correlations and the common offset of host phase - ZT (`Dataset.timed`); tissue offsets arg z of the
anchored fit (`Dataset.fit`, all tissues anchors, leave one out) against -arg g of the per-tissue timed cosinor gain (14 genes
with median CPM >= 5, regressed on the global GTEx waveform). Leave-one-animal-out host (stability). Sensitivity: peripheral and
neural anchors alone. Null: ZT (or tissue phase) permuted, 10000 permutations (--subsample: 500).
Human gains: fast deaths (Hardy 1-2), cosinor with covariates (sex, age bracket, Hardy, ischaemic time and its square) per
matched GTEx tissue, local aggregate; no per-donor values are written.

    python compute_baboon_truth.py [--subsample N]   (N: random animals; the 12 evenly spaced ZTs are the design)
Outputs: paper/results/figure_1/output/method/baboon_{readouts,animals,tissues,pairs,stability}.csv, baboon_gate.csv
Legacy: legacy/desync/validation/baboon/scripts/{prep,score,human_gain}.py (output/readouts.csv, tissues.csv, animals.csv)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402

FOLDER = "method"
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
NEURAL = {"AMY", "ARC", "CER", "DMH", "HAB", "HIP", "LGP", "LH", "MGP", "MMB", "OLB", "PON", "PRA", "PRC", "PUT", "PVN", "SCN",
          "SON", "SUN", "THA", "VIC", "VMH", "PIN", "PIT", "RET", "RPE", "COR", "ONH", "IRI"}


def timed_gain(sub: pd.DataFrame) -> tuple[complex, int]:
    """local: not in clocksync. Per-tissue 12-ZT cosinor of log(CPM + 0.5) on the clock genes with median CPM >= 5, regressed on the
    waveform: complex gain g."""
    cpm = sub[S.CLOCK].to_numpy(float) / sub.depth.to_numpy(float)[:, None] * 1e6
    ex = np.median(cpm, 0) >= 5
    Y = np.log(cpm[:, ex] + 0.5)
    t = sub.zt.to_numpy(float) * W
    D = np.c_[np.ones(len(t)), np.cos(t), np.sin(t)]
    B = np.linalg.lstsq(D, Y, rcond=None)[0]
    c = B[1] + 1j * B[2]
    x = (S.WF.amp * np.exp(1j * S.WF.psi))[ex]
    return complex((np.conj(x) @ c) / (np.conj(x) @ x)), int(ex.sum())


def human_gains(tissues) -> pd.DataFrame:
    """local: not in clocksync. GTEx timed tissue gain (fast deaths, covariates sex, age, Hardy, ischaemic time and square) with a
    donor-bootstrap SE of its phase. Aggregates only."""
    rows = []
    for tis in tissues:
        pb, d = gtex.load(tissues=[tis], hardy=(1, 2), depth="library_size")
        pb = pb.drop_duplicates("donor").set_index("donor")
        d = d.loc[pb.index]
        cpm = pb[S.CLOCK].to_numpy(float) / pb.depth.to_numpy(float)[:, None] * 1e6
        t = d.tod.to_numpy(float) % 24
        isch = d.ischaemic_h.fillna(d.ischaemic_h.median()).to_numpy(float)
        isch = (isch - isch.mean()) / (isch.std() + 1e-9)
        Z = np.c_[pd.get_dummies(d.age.astype(str), drop_first=True).to_numpy(float), d.male.to_numpy(float),
                  pd.get_dummies(d.hardy.astype(str), drop_first=True).to_numpy(float), isch, isch ** 2]
        x_all = S.WF.amp * np.exp(1j * S.WF.psi)

        def gain(i):
            ex = np.median(cpm[i], 0) >= 5
            Y = np.log(cpm[i][:, ex] + 0.5)
            D = np.c_[np.ones(len(i)), Z[i], np.cos(W * t[i]), np.sin(W * t[i])]
            B = np.linalg.lstsq(D, Y, rcond=None)[0]
            c = B[-2] + 1j * B[-1]
            xc = x_all[ex]
            return (np.conj(xc) @ c) / (np.conj(xc) @ xc), int(ex.sum())

        g, ng = gain(np.arange(len(t)))
        rng = np.random.default_rng(0)
        bs = []
        for _ in range(50):
            try:
                bs.append(gain(rng.choice(len(t), len(t)))[0])
            except np.linalg.LinAlgError:
                continue
        dphi = np.angle(np.asarray(bs) * np.conj(g)) / W
        rows.append(dict(tissue=tis, n=len(t), clock_genes=ng, abs_g=abs(g), phi_g_h=(np.angle(g) / W) % 24,
                         phi_g_se_h=float(np.std(dphi, ddof=1))))
    return pd.DataFrame(rows)


def main():
    ap = cli.parser(__doc__)
    args = ap.parse_args()
    opt = cli.setup(args)
    out = out_dir(FOLDER, opt["quick"])
    n_perm = 500 if opt["quick"] else 10000
    pb, donors = data._builder("baboon").build()
    donors.index = donors.index.astype(str)
    pb["donor"] = pb.donor.astype(str)
    if args.subsample is not None:
        keep = np.random.default_rng(args.seed).permutation(sorted(pb.donor.unique()))[:args.subsample]
        pb, donors = pb[pb.donor.isin(keep)].reset_index(drop=True), donors.loc[donors.index.isin(keep)]
    zt = donors.zt.astype(float)
    kept = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb, anchors=kept, donors=donors, name="baboon")
    print(f"A36: {len(kept)} tissues x {pb.donor.nunique()} animals", flush=True)

    # readout 1: host phase vs ZT (+ sensitivity: peripheral / neural anchors)
    rows = []
    r1 = ds.timed(zt, n_perm=n_perm)
    rows.append(dict(readout="1 PRIMARY host phase vs ZT (all tissues)", **r1.to_dict()))
    for lab, anc in (("peripheral anchors", [k for k in kept if k not in NEURAL]),
                     ("neural/ocular anchors", [k for k in kept if k in NEURAL])):
        if len(anc) >= 2:
            rows.append(dict(readout=f"1 sensitivity: {lab} ({len(anc)})", **ds.timed(zt, n_perm=n_perm, anchors=anc).to_dict()))
    # animals and stability: leave one animal out
    fit = ds.fit(n_boot=0, gate=False)
    host = fit.host
    phi = np.angle(host.x + 1j * host.y)
    animals = pd.DataFrame(dict(animal=host.index, zt=zt.reindex(host.index).to_numpy(), host_phase_h=(phi / W) % 24,
                                host_amp=np.hypot(host.x, host.y)))
    animals["err_h"] = S.cdiff(animals.host_phase_h - r1.offset_h, animals.zt)
    animals.to_csv(out / "baboon_animals.csv", index=False)
    lo = []
    for a in host.index:
        sub = ds.subset(donors=[x for x in host.index if x != a])
        h2 = sub.fit(n_boot=0, gate=False).host
        r = sub.timed(zt, n_perm=200)
        common = h2.index
        sim = abs(cs.stats.ccorr(np.c_[h2.x, h2.y], np.c_[host.loc[common, "x"], host.loc[common, "y"]]))
        lo.append(dict(left_out=a, r_js=r.r_js, r_fl=r.r_fl, R=r.R, sim_to_full=sim))
    lo = pd.DataFrame(lo)
    lo.to_csv(out / "baboon_stability.csv", index=False)
    print(f"LOAO Fisher-Lee min {lo.r_fl.min():.3f}; R min {lo.R.min():.3f}; similarity min {lo.sim_to_full.min():.3f}", flush=True)

    # readout 2: untimed tissue offsets vs timed tissue phases
    T = fit.table[["coupling", "log_amp_rel", "offset_h", "boundary"]].copy()
    d = ds.data
    gains, bph = {}, {}
    for k in kept:
        sub = pb[pb.compartment == k].assign(zt=lambda x: zt.reindex(x.donor).to_numpy())
        g, ng = timed_gain(sub)
        gains[k] = (g, ng)
        i = d.compartment == k
        bph[k] = S.timed_loading(d.beta[i], zt.reindex(d.resection[i]).to_numpy())
    T["neural"] = [k in NEURAL for k in T.index]
    T["g_abs"] = [abs(gains[k][0]) for k in T.index]
    T["g_arg"] = [np.angle(gains[k][0]) for k in T.index]
    T["g_genes"] = [gains[k][1] for k in T.index]
    T["beta_amp"] = [abs(bph[k]) for k in T.index]
    T["beta_arg"] = [np.angle(bph[k]) for k in T.index]
    T["timed_phase"] = -T.g_arg
    tc = np.angle(np.exp(1j * T.timed_phase).mean())
    T["timed_offset_h"] = np.angle(np.exp(1j * (T.timed_phase - tc))) / W
    T["g_phi_h"] = (T.g_arg / W) % 24
    T["z_arg"] = T.offset_h * W
    T.to_csv(out / "baboon_tissues.csv")
    ag = lambda a, b: cs.timed.agreement(a, b, n_perm=n_perm)
    rows.append(dict(readout="2 untimed offset arg z vs timed -arg g (all tissues)", **ag(T.z_arg.to_numpy(), T.timed_phase.to_numpy()).to_dict()))
    for lab, m in (("peripheral", ~T.neural), ("neural/ocular", T.neural)):
        if m.sum() >= 4:
            rows.append(dict(readout=f"2 sensitivity: {lab}", **ag(T.z_arg[m].to_numpy(), T.timed_phase[m].to_numpy()).to_dict()))
    rows.append(dict(readout="2 check (not a readout): -arg g vs timed beta phase", **ag(T.timed_phase.to_numpy(), T.beta_arg.to_numpy()).to_dict()))

    # readout 3: cross-species
    pairs = {b: g for b, g in PAIRS.items() if b in T.index}
    hg = human_gains(sorted(set(pairs.values()))).set_index("tissue")
    P = pd.DataFrame([dict(baboon=b, gtex=g) for b, g in pairs.items() if g in hg.index])
    P["baboon_g_phi_h"] = T.g_phi_h.reindex(P.baboon).to_numpy()
    P["baboon_g_abs"] = T.g_abs.reindex(P.baboon).to_numpy()
    P["baboon_untimed_lag_h"] = ((-T.z_arg / W) % 24).reindex(P.baboon).to_numpy()
    P["human_g_phi_h"] = hg.phi_g_h.reindex(P.gtex).to_numpy()
    P["human_g_abs"] = hg.abs_g.reindex(P.gtex).to_numpy()
    P["human_phi_se_h"] = hg.phi_g_se_h.reindex(P.gtex).to_numpy()
    P["human_n"] = hg.n.reindex(P.gtex).to_numpy()
    P.to_csv(out / "baboon_pairs.csv", index=False)
    if len(P) >= 6:
        bw, hw = P.baboon_g_phi_h.to_numpy() * W, P.human_g_phi_h.to_numpy() * W
        per = ~P.baboon.isin(BRAIN_PAIRS).to_numpy()
        rows.append(dict(readout="3 cross-species: baboon timed arg g vs human arg g", **ag(bw, hw).to_dict()))
        rows.append(dict(readout="3 sensitivity: peripheral pairs", **ag(bw[per], hw[per]).to_dict()))
        ul = P.baboon_untimed_lag_h.to_numpy() * W
        rows.append(dict(readout="3 sensitivity: baboon untimed (-arg z) vs human arg g", **ag(ul, hw).to_dict()))
        rows.append(dict(readout="3 sensitivity: baboon untimed, peripheral pairs", **ag(ul[per], hw[per]).to_dict()))
    R = pd.DataFrame(rows)
    R.to_csv(out / "baboon_readouts.csv", index=False)
    pd.set_option("display.width", 250)
    print(R.round(4).to_string(index=False))
    p1 = R.iloc[0]
    print(f"PRIMARY as written (JS r > 0.5, p < 0.05): {'PASSES' if p1.r_js > 0.5 and p1.p_js < 0.05 else 'FAILS'} (r {p1.r_js:.3f}, p {p1.p_js:.4g})")
    print(f"PRIMARY with Fisher-Lee (post hoc): {'PASSES' if p1.r_fl > 0.5 and p1.p_fl < 0.05 else 'FAILS'} (r_fl {p1.r_fl:.3f}, p {p1.p_fl:.4g})")
    gt = ds.gate(n_perm=50 if opt["quick"] else 200)
    gt.to_series().to_csv(out / "baboon_gate.csv", header=False)
    print(gt.summary())


if __name__ == "__main__":
    main()
