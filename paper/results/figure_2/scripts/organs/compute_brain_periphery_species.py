"""A42 (+ A103.3): is the brain partly separate from the periphery outside humans? Sheep (untimed) and baboon (timed).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A42; A103 for H103.3). Tissue groups fixed beforehand: brain = central
nervous system regions only; pineal, pituitary and ocular tissues excluded from both groups in the primary (sensitivity:
counted as brain); periphery = all other tissues.
  Readout 1 (message 2): two untimed free hosts on the same animals, one from the brain tissues only and one from the
    peripheral tissues only; r* = |complex correlation| of the two hosts over animals, disattenuated by each host's
    split-half reliability (`Dataset.agreement`, 50 splits); animal bootstrap interval. Pass (partly separate): r* < 1 and
    the upper bootstrap bound < 1. Sheep needs >= 20 animals with >= 2 brain and >= 2 peripheral tissues; the study
    (BioProject) overlap of the two blocks is reported, with sensitivities on animals whose brain and peripheral samples
    come from identical study sets. Baboon (12 animals, ZT known): also the timed version, descriptive (host phase minus
    ZT in each block; unit-circle, Fisher-Lee and rotated-state correlations, 10,000 permutations).
  Readout 2 (message 2): coupling of each tissue to the peripheral-anchor host (leave one out for peripheral tissues); brain
    vs periphery by one-sided Mann-Whitney. Pass: p < 0.05 in each species.
  Readout 3 (message 3): tissue couplings on the all-tissue host against GTEx fast-death couplings of the matched human
    tissues (Spearman, one-sided permutation p over tissues, 10,000). Sheep: primary, >= 10 matched tissues, rho > 0 and
    p < 0.05; baboon: secondary.
  A103.3: readout 1 in sheep with the clock plane without PER1 and PER2 (12-gene waveform).

    python compute_brain_periphery_species.py [--subsample N] [--jobs 1] [--n-boot 200] [--n-boot-sens 100]
Outputs: paper/results/figure_2/output/organs/{r1_sheep_readout,r1_study_overlap,r1_baboon_readout,r1_baboon_timed,r2_couplings,
r2_readout,r3_readout,r3_pairs,a103_sheep_readout}.csv
Legacy: legacy/desync/validation/cross_species/scripts/{r1,r23,common,human_ext}.py,
legacy/desync/validation/brain_artefact/scripts/a103_sheep.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402

FOLDER = "organs"
W = 2 * np.pi / 24
N_PERM = 10_000
SHEEP_BRAIN = ["Cerebellum", "Hypothalamus", "Cerebrum", "Brainstem", "Medulla oblongata", "Cerebral cortex",
               "Hippocampus", "Optic chiasm"]
SHEEP_EXCL = ["Pineal", "Pituitary"]
BAB_BRAIN = ["AMY", "ARC", "CER", "DMH", "HAB", "HIP", "LGP", "LH", "MGP", "MMB", "OLB", "PON", "PRA", "PRC", "PUT",
             "PVN", "SCN", "SON", "SUN", "THA", "VIC", "VMH"]
BAB_EXCL = ["PIN", "PIT", "RET", "RPE", "COR", "ONH"]
# sheep compartment -> GTEx tissue (A34, fixed before scoring) and the secondary additions
MATCH = {
    "Liver": "Liver", "Muscle": "Muscle - Skeletal", "Kidney": "Kidney - Cortex", "Cerebellum": "Brain - Cerebellum",
    "Hypothalamus": "Brain - Hypothalamus", "Lung": "Lung", "Colon": "Colon - Transverse",
    "Heart": "Heart - Left Ventricle", "Pituitary": "Pituitary", "Subcutaneous adipose": "Adipose - Subcutaneous",
    "Perirenal adipose": "Adipose - Visceral (Omentum)", "Thyroid": "Thyroid", "Artery": "Artery - Aorta",
    "Skin": "Skin - Not Sun Exposed (Suprapubic)", "Ovary": "Ovary", "Uterus": "Uterus",
    "Cerebral cortex": "Brain - Cortex", "Hippocampus": "Brain - Hippocampus", "Adrenal": "Adrenal Gland",
    "Mammary gland": "Breast - Mammary Tissue", "Abomasum": "Stomach",
}
SECONDARY = {"Rumen": "Esophagus - Mucosa", "Reticulum": "Esophagus - Mucosa", "Omasum": "Esophagus - Mucosa",
             "Cerebrum": "Brain - Frontal Cortex (BA9)"}
# baboon tissue -> GTEx tissue (A36)
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


def groups(species: str, comps, sens: bool = False):
    """(brain, periphery, excluded) tissue lists present in `comps` (pineal, pituitary, ocular excluded in the primary)."""
    b, x = (SHEEP_BRAIN, SHEEP_EXCL) if species == "sheep" else (BAB_BRAIN, BAB_EXCL)
    comps = list(comps)
    assert not (set(b) | set(x)) - set(comps), (set(b) | set(x)) - set(comps)
    brain = [c for c in comps if c in b or (sens and c in x)]
    excl = [] if sens else [c for c in comps if c in x]
    per = [c for c in comps if c not in b and c not in x]
    return brain, per, excl


def eligible(table: pd.DataFrame, brain, per):
    g = pd.DataFrame(dict(a=table.resection, b=table.compartment.isin(brain), p=table.compartment.isin(per))).groupby("a")[["b", "p"]].sum()
    return list(g.index[(g.b >= 2) & (g.p >= 2)])


def readout1(label, ds, brain, per, n_boot, n_split, min_common, rows, boots=None):
    """`Dataset.agreement` of the brain and the peripheral host, as a row of the readout table."""
    r = ds.agreement({"brain": brain, "periphery": per}, n_boot=n_boot, n_split=n_split, seed=0, min_common=min_common)
    row = dict(analysis=label, n_brain_tissues=len(brain), n_periph_tissues=len(per), **r.rename(
        {"rel_brain": "rel_brain", "rel_periphery": "rel_periph", "n_donors": "n_animals", "n_eligible": "n_eligible"}).to_dict())
    row["pass"] = bool(r.get("separate", False)) if n_boot else np.nan
    rows.append(row)
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    return row


def dataset(pb, genes=None, anchors=None, mixed=False, name=""):
    wf = cs.gtex() if genes is None else cs.gtex().subset(genes)
    cols = [c for c in ("donor", "compartment", "depth", "stratum") if c in pb] + list(wf.genes)
    return cs.Dataset(pb[cols], anchors=sorted(pb.compartment.unique()) if anchors is None else anchors, waveform=wf,
                      mixed_strata=mixed, name=name, min_donors=5)


def sheep_r1(pb, args, run, out):
    ds = dataset(pb, mixed=True, name="sheep")
    comps = sorted(ds.compartments)
    rows = []
    n_split = 50 if not run["quick"] else 6
    nb, nbs = (args.n_boot, args.n_boot_sens) if not run["quick"] else (3, 2)
    brain, per, _ = groups("sheep", comps)
    el = eligible(ds.table, brain, per)
    t = ds.table[ds.table.resection.isin(el) & ds.table.compartment.isin(brain + per)]
    print(f"PRIMARY: {len(el)} eligible animals", flush=True)
    if len(el) < 20 and not run["quick"]:
        rows.append(dict(analysis="PRIMARY", status=f"infeasible ({len(el)} < 20 animals)"))
    else:
        # study overlap of the two blocks (reported: a brain-only study effect could pose as separation)
        tt = pd.DataFrame(dict(a=t.resection, s=t.stratum, b=t.compartment.isin(brain)))
        ov = tt.groupby(["s", "b"]).a.nunique().unstack(fill_value=0).rename(columns={True: "brain_animals",
                                                                                    False: "periph_animals"})
        sb, sp = tt[tt.b].groupby("a").s.apply(set), tt[~tt.b].groupby("a").s.apply(set)
        ov["animals_with_both_in_this_study"] = [sum((s in sb[a]) and (s in sp[a]) for a in el) for s in ov.index]
        ov.to_csv(out / "r1_study_overlap.csv")
        same = [a for a in el if sb[a] == sp[a]]
        print(ov.to_string(), f"\nanimals sharing >= 1 study: {sum(len(sb[a] & sp[a]) > 0 for a in el)}/{len(el)}; "
              f"identical study sets: {len(same)}", flush=True)
        readout1("PRIMARY", ds, brain, per, nb, n_split, 10, rows)
        if len(same) >= 12:
            readout1(f"sens: identical brain/periphery study sets ({len(same)} animals)", ds.subset(donors=same), brain, per,
                     nbs, n_split, 10, rows)
    b2, p2, _ = groups("sheep", comps, True)
    readout1("sens: pineal/pituitary as brain", ds, b2, p2, nbs, n_split, 10, rows)
    pd.DataFrame(rows).to_csv(out / "r1_sheep_readout.csv", index=False)
    return ds


def baboon_timed(ds, brain, per, el, label, rng, n_perm):
    """Descriptive timed version: departure of each block's host from ZT (host phase minus 2 pi ZT / 24), per animal."""
    t = ds.table[ds.table.resection.isin(el)]

    def block(comps):
        sub = t[t.compartment.isin(comps)].rename(columns={"resection": "donor"}).reset_index(drop=True)
        d = cs.Dataset(sub, anchors=comps, min_donors=5)
        h = d.fit(n_boot=0, gate=False, targets=[comps[0]]).host
        f = h[["x", "y"]].to_numpy() - h[["x", "y"]].to_numpy().mean(0)
        zt = np.array([int(a[2:]) for a in h.index], float) * W
        z = (f[:, 0] + 1j * f[:, 1]) * np.exp(-1j * zt)
        return pd.Series(np.angle(z), index=h.index), pd.Series(z, index=h.index)

    db, wb = block(brain)
    dp, wp = block(per)
    dp, wp = dp[db.index], wp[db.index]
    perms = [rng.permutation(len(dp)) for _ in range(n_perm)]
    cz = lambda a, b: abs(cs.stats.ccorr(np.c_[a.real, a.imag], np.c_[b.real, b.imag]))
    u = lambda x: np.c_[np.cos(x), np.sin(x)]
    r_u = abs(cs.stats.ccorr(u(db), u(dp)))
    r_fl = cs.stats.fisher_lee(db.to_numpy(), dp.to_numpy())
    nu = [abs(cs.stats.ccorr(u(db), u(dp.to_numpy()[p]))) for p in perms]
    nfl = [cs.stats.fisher_lee(db.to_numpy(), dp.to_numpy()[p]) for p in perms]
    nw = [cz(wb.to_numpy(), wp.to_numpy()[p]) for p in perms]
    r_w = cz(wb.to_numpy(), wp.to_numpy())
    pv = lambda null, o: cs.stats.perm_p(null, o, "greater", tie=0.0)
    csd = lambda x: np.sqrt(-2 * np.log(abs(np.exp(1j * x).mean()))) / W
    return dict(analysis=label, n=len(db), r_unit_departure=r_u, p_unit=pv(nu, r_u), r_fisher_lee=r_fl, p_fl=pv(nfl, r_fl),
                r_rotated_state=r_w, p_rotated=pv(nw, r_w), R_brain=abs(np.exp(1j * db).mean()),
                R_periph=abs(np.exp(1j * dp).mean()), circ_sd_brain_h=csd(db.to_numpy()), circ_sd_periph_h=csd(dp.to_numpy()))


def baboon_r1(pb, args, run, out):
    ds = dataset(pb, name="baboon")
    comps = sorted(ds.compartments)
    rows, trows = [], []
    rng = np.random.default_rng(0)
    nb, nbs = (args.n_boot, args.n_boot_sens) if not run["quick"] else (3, 2)
    n_split = 50 if not run["quick"] else 6
    n_perm = N_PERM if not run["quick"] else 500
    for sens in (False, True):
        brain, per, _ = groups("baboon", comps, sens)
        el = eligible(ds.table, brain, per)
        lab = "sens: pineal/pituitary/ocular as brain" if sens else "PRIMARY"
        print(f"{lab}: {len(el)} eligible animals, {len(brain)} brain, {len(per)} periphery", flush=True)
        readout1(lab, ds, brain, per, nbs if sens else nb, n_split, len(el), rows)
        trows.append(baboon_timed(ds, brain, per, el, lab, rng, n_perm))
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in trows[-1].items()}, flush=True)
    pd.DataFrame(rows).to_csv(out / "r1_baboon_readout.csv", index=False)
    pd.DataFrame(trows).to_csv(out / "r1_baboon_timed.csv", index=False)
    return ds


def readout2(species, pb, comps, mixed):
    brain, per, excl = groups(species, comps)
    ds = dataset(pb, anchors=per, mixed=mixed, name=species)
    tab = ds.fit(n_boot=0, gate=False, targets=brain + per + excl).table
    tab["group"] = ["brain" if k in brain else "periphery" if k in per else "excluded" for k in tab.index]
    tab.insert(0, "species", species)
    rows = []
    for lab, b in (("primary", brain), ("sensitivity: pineal/pituitary/ocular as brain", brain + excl)):
        cb, cp = tab.coupling.reindex(b).dropna().to_numpy(), tab.coupling.reindex(per).dropna().to_numpy()
        mw = stats.mannwhitneyu(cb, cp, alternative="less")
        rows.append(dict(species=species, analysis=lab, n_brain=len(cb), n_periph=len(cp), median_brain=np.median(cb),
                         median_periph=np.median(cp), U=mw.statistic, p=mw.pvalue, **{"pass": bool(mw.pvalue < 0.05)}))
    return tab, rows


def spearman_perm(x, y, n_perm, seed=0):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rho = stats.spearmanr(x, y).statistic
    rng = np.random.default_rng(seed)
    null = np.array([stats.spearmanr(x, rng.permutation(y)).statistic for _ in range(n_perm)])
    return rho, cs.stats.perm_p(null, rho, "greater", tie=0.0)


def r3(species, analysis, own, human, pairs, pair_rows, crit, n_perm):
    m = [(k, v) for k, v in pairs if np.isfinite(own.get(k, np.nan)) and np.isfinite(human.get(v, np.nan))]
    x, y = [own[k] for k, _ in m], [human[v] for _, v in m]
    rho, p = spearman_perm(x, y, n_perm) if len(m) >= 3 else (np.nan, np.nan)
    for k, v in m:
        pair_rows.append(dict(species=species, analysis=analysis, tissue=k, human=v, coupling=own[k], human_coupling=human[v]))
    return {"species": species, "analysis": analysis, "n": len(m), "rho": rho, "p": p,
            "pass": bool(len(m) >= 10 and rho > 0 and p < 0.05) if crit else np.nan}


def human_control(tissues):
    """Anchored couplings of GTEx fast-death tissues (A34 `human_truth.control`): every tissue an anchor."""
    pb, _ = gtex.load(tissues=list(tissues), hardy=gtex.FAST, genes="clock", depth="library_size", require_tod=False)
    ds = dataset(pb.assign(donor=pb.donor), anchors=sorted(pb.compartment.unique()))
    return ds.fit(n_boot=0, gate=False).table.coupling


def human_tissues(extra=()):
    """Tissues of the A34 control set (matched tissues and the 15 reference tissues with >= 20 fast-death donors that
    have a time of death)."""
    cand = sorted(set(MATCH.values()) | set(S.V2_TISSUES) | set(extra))
    pb, _ = gtex.load(tissues=cand, hardy=gtex.FAST, genes="clock", depth="library_size", require_tod=True)
    n = pb.groupby("compartment").donor.nunique()
    return sorted(n.index[n >= 20])


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--n-boot-sens", type=int, default=100)
    ap.add_argument("--skip-r1", action="store_true", help="only readouts 2 and 3 (cheap; r1 and A103.3 are the slow part)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = N_PERM if not run["quick"] else 500
    pd.set_option("display.width", 220)

    spb, _ = data._builder("sheep").build()
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(spb.donor.unique()))[:args.subsample]
        spb = spb[spb.donor.isin(keep)]
    bpb, _ = data._builder("baboon").build()
    sheep_comps, bab_comps = sorted(spb.compartment.unique()), sorted(bpb.compartment.unique())

    # ---- readout 1 and A103.3
    if not args.skip_r1:
        sheep_r1(spb, args, run, out)
        baboon_r1(bpb, args, run, out)
    g12 = [g for g in S.CLOCK if g not in ("PER1", "PER2")]
    brain, per, _ = groups("sheep", sheep_comps)
    rows = []
    if not args.skip_r1:
        ds12 = dataset(spb, genes=g12, mixed=True, name="sheep12")
        readout1("PRIMARY H103.3: no PER1/PER2", ds12, brain, per, args.n_boot if not run["quick"] else 3,
                 50 if not run["quick"] else 6, 10, rows)
        pd.DataFrame(rows).to_csv(out / "a103_sheep_readout.csv", index=False)

    # ---- readouts 2 and 3
    C2, R2, R3, P3 = [], [], [], []
    c2, r2 = readout2("sheep", spb, sheep_comps, True)
    C2.append(c2)
    R2 += r2
    c2, r2 = readout2("baboon", bpb, bab_comps, False)
    C2.append(c2)
    R2 += r2
    own_sheep = dataset(spb, mixed=True).fit(n_boot=0, gate=False).table.coupling.to_dict()
    own_bab = dataset(bpb).fit(n_boot=0, gate=False).table.coupling.to_dict()
    hA = human_control(human_tissues()).to_dict()
    hX = human_control(human_tissues(extra=PAIRS.values())).to_dict()
    brain_s, per_s, excl_s = groups("sheep", sheep_comps)
    nonbrain = lambda k, v: k not in brain_s + excl_s and not v.startswith("Brain") and v != "Pituitary"
    hAX = hA | {k: v for k, v in hX.items() if k not in hA}
    R3.append(r3("sheep", "PRIMARY: MATCH, A34 human control", own_sheep, hA, list(MATCH.items()), P3, True, n_perm))
    R3.append(r3("sheep", "sens: MATCH + SECONDARY", own_sheep, hAX, list((MATCH | SECONDARY).items()), P3, False, n_perm))
    R3.append(r3("sheep", "sens: non-brain matched only", own_sheep, hA,
                 [(k, v) for k, v in MATCH.items() if nonbrain(k, v)], P3, False, n_perm))
    R3.append(r3("sheep", "sens: MATCH, extended human fit", own_sheep, hX, list(MATCH.items()), P3, False, n_perm))
    R3.append(r3("baboon", "secondary: PAIRS, extended human fit", own_bab, hX, list(PAIRS.items()), P3, True, n_perm))
    R3.append(r3("baboon", "sens: peripheral pairs", own_bab, hX,
                 [(k, v) for k, v in PAIRS.items() if k not in BRAIN_PAIRS and k != "PIT"], P3, False, n_perm))
    C2 = pd.concat(C2)
    C2.to_csv(out / "r2_couplings.csv")
    pd.DataFrame(R2).to_csv(out / "r2_readout.csv", index=False)
    pd.DataFrame(R3).to_csv(out / "r3_readout.csv", index=False)
    pd.DataFrame(P3).to_csv(out / "r3_pairs.csv", index=False)
    print(pd.DataFrame(R2).round(4).to_string(index=False), "\n", pd.DataFrame(R3).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
