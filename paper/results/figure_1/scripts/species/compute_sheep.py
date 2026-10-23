"""A34: farm-animal GTEx bulk. Is the order of tissue clock phases conserved between sheep and human?

Pre-registered (A34, written before running; its run note records the feasibility result). SheepGTEx v0 is the only farm species
with animal IDs (343 animals, 45 tissues / cell types, 38 BioProjects as strata; no slaughter times, cattle and pig infeasible).
Readouts: (a) the clock gate over all tissues (and per tissue: the share of the tissue's pairs; per BioProject), (b, c) the
anchored host with all tissues as anchors, tissue offsets, amplitudes and couplings with an animal bootstrap, (e) primary:
circular correlation (Jammalamadaka-SenGupta, one-sided permutation p, 10000) of the sheep tissue offsets with the human GTEx timed
tissue phase d_t = -arg g_t over the matched tissues (>= 6; r > 0 with p < 0.05), jackknife range, bootstrap interval of r.
Secondary: forestomach and cerebrum matches. Positive control: the same untimed-vs-timed comparison on GTEx fast deaths themselves
(tissues as anchors). Exploratory (after the primary): Spearman of the sheep tissue amplitude with the human timed gain.

    python compute_sheep.py [--subsample N] [--n-boot B] [--jobs J] [--tissue-gates]    (N: animals)
Outputs: paper/results/figure_1/output/species/a34_{gate,offsets,primary,human_gain,human_control,amplitude}.csv
Legacy: validation/farm_gtex/scripts/{fit_sheep,human_truth,match}.py (per-tissue and per-study gate re-implementations: `--tissue-gates`
uses `Dataset.gate(targets=, hosts=)`)
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
from common import cli, data, gtex  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "species"
W, HOURS = S.W, S.HOURS
MATCH = {
    "Liver": "Liver", "Muscle": "Muscle - Skeletal", "Kidney": "Kidney - Cortex", "Cerebellum": "Brain - Cerebellum",
    "Hypothalamus": "Brain - Hypothalamus", "Lung": "Lung", "Colon": "Colon - Transverse", "Heart": "Heart - Left Ventricle",
    "Pituitary": "Pituitary", "Subcutaneous adipose": "Adipose - Subcutaneous", "Perirenal adipose": "Adipose - Visceral (Omentum)",
    "Thyroid": "Thyroid", "Artery": "Artery - Aorta", "Skin": "Skin - Not Sun Exposed (Suprapubic)", "Ovary": "Ovary", "Uterus": "Uterus",
    "Cerebral cortex": "Brain - Cortex", "Hippocampus": "Brain - Hippocampus", "Adrenal": "Adrenal Gland",
    "Mammary gland": "Breast - Mammary Tissue", "Abomasum": "Stomach",
}
SECONDARY = {"Rumen": "Esophagus - Mucosa", "Reticulum": "Esophagus - Mucosa", "Omasum": "Esophagus - Mucosa",
             "Cerebrum": "Brain - Frontal Cortex (BA9)"}
CONTROL_EXTRA = ["Artery - Tibial", "Esophagus - Mucosa", "Heart - Atrial Appendage", "Nerve - Tibial",
                 "Skin - Sun Exposed (Lower leg)", "Whole Blood"]


def circ_corr(a, b):
    return float(cs.stats.js_corr(a, b))


def primary(off, human, match, n_perm, seed=0):
    m = [(k, v) for k, v in match.items() if k in off.index and v in human.index and np.isfinite(off.offset_h[k])]
    if len(m) < 6:
        return dict(n=len(m), r=np.nan, p=np.nan, status="infeasible (< 6 matched tissues)"), m
    a = np.array([off.offset_h[k] for k, _ in m]) / HOURS
    b = np.array([human.d_h[v] for _, v in m]) / HOURS
    r = circ_corr(a, b)
    rng = np.random.default_rng(seed)
    null = np.array([circ_corr(a, rng.permutation(b)) for _ in range(n_perm)])
    jk = [circ_corr(np.delete(a, i), np.delete(b, i)) for i in range(len(a))]
    return dict(n=len(m), r=r, p=(1 + np.sum(null >= r)) / (1 + n_perm), null_q95=float(np.quantile(null, .95)),
                jackknife_min=float(np.min(jk)), jackknife_max=float(np.max(jk)), status="ok"), m


def _boot(ds, seed):
    limit_threads(1)
    try:
        return ds.resample(np.random.default_rng(seed)).fit(n_boot=0, gate=False).table.offset_h
    except (ValueError, np.linalg.LinAlgError):
        return None


def tissue_gate(ds, k, n_perm):
    others = [c for c in ds.compartments if c != k]
    g = ds.gate(n_perm=n_perm, targets=[k], hosts=others)
    return dict(compartment=k, clock_share=g.clock_share, null_mean=g.share_null_mean, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase,
                n_pairs=g.n_pairs)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=None)
    ap.add_argument("--tissue-gates", action="store_true")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    out = out_dir(FOLDER, quick)
    n_boot = args.n_boot if args.n_boot is not None else (5 if quick else 100)
    n_perm = 1000 if quick else 10_000
    pb, donors = data.load("sheep", subsample=args.subsample, seed=args.seed)
    pb = pb.copy()
    ds = cs.Dataset(pb, anchors=None, donors=donors, name="sheep", mixed_strata=True)
    comps = ds.compartments
    print(f"{pb.donor.nunique()} animals, {len(comps)} compartments, {pb.stratum.nunique()} strata", flush=True)
    g = ds.gate(n_perm=200 if quick else 1000)
    pd.DataFrame([g.to_series()]).to_csv(out / "a34_gate.csv", index=False)
    print(g.summary())
    if args.tissue_gates and not quick:
        pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(tissue_gate)(ds, k, 200) for k in comps)).to_csv(out / "a34_gate_tissue.csv", index=False)
    f = ds.fit(n_boot=0, gate=False)
    tab = f.table.copy()
    reps = [r for r in Parallel(n_jobs=args.jobs)(delayed(_boot)(ds, s) for s in range(n_boot)) if r is not None]
    if reps:
        B = pd.DataFrame(reps)
        tab["offset_sd_h"] = [np.std(np.angle(np.exp(1j * (B[k].dropna() - tab.offset_h[k]) / HOURS)) * HOURS) for k in tab.index]
    tab["n_boot"] = len(reps)
    tab.to_csv(out / "a34_offsets.csv")

    # human truths: timed GTEx gains of the matched tissues (fast deaths) -> d_t = -arg g_t (h)
    tis = sorted(set(MATCH.values()) | set(SECONDARY.values()) | set(CONTROL_EXTRA))
    hg = S.gtex_gains(tis, n_boot=0 if quick else 200, subsample=args.subsample, seed=args.seed)
    hg["d_h"] = (-np.angle(hg.g_re + 1j * hg.g_im) * HOURS)
    hg.to_csv(out / "a34_human_gain.csv", index=False)
    hum = hg.set_index("tissue")
    rows = []
    for lab, match in (("primary", MATCH), ("secondary: + forestomach, cerebrum", MATCH | SECONDARY)):
        r, m = primary(tab, hum, match, n_perm)
        r.update(analysis=lab, tissues=";".join(k for k, _ in m))
        if lab == "primary" and r["status"] == "ok" and reps:  # bootstrap interval of r (offsets resampled)
            rb = []
            for rep in reps:
                x = pd.DataFrame(dict(offset_h=rep))
                a = np.array([x.offset_h.get(k, np.nan) for k, _ in m]) / HOURS
                b = np.array([hum.d_h[v] for _, v in m]) / HOURS
                ok = np.isfinite(a)
                if ok.sum() >= 6:
                    rb.append(circ_corr(a[ok], b[ok]))
            r["r_boot_lo"], r["r_boot_hi"] = (np.quantile(rb, .025), np.quantile(rb, .975)) if rb else (np.nan, np.nan)
        rows.append(r)
    # positive control: the same comparison inside GTEx (fast deaths; tissues as anchors)
    ctrl_t = sorted(set(MATCH.values()) | set(CONTROL_EXTRA))
    pbg, dg = gtex.load(ctrl_t, depth="library_size")
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(pbg.donor.unique()))[:max(args.subsample, 60)]
        pbg = pbg[pbg.donor.isin(keep)]
    dsg = cs.Dataset(pbg, anchors=None, name="gtex_fast")
    tg = dsg.fit(n_boot=0, gate=False).table
    tg.to_csv(out / "a34_human_control_offsets.csv")
    for lab, tl in (("human control: matched tissues", sorted(set(MATCH.values()))), ("human control: all fitted tissues", ctrl_t)):
        tl = [t for t in tl if t in tg.index and t in hum.index]
        a = tg.offset_h[tl].to_numpy() / HOURS
        b = hum.d_h[tl].to_numpy() / HOURS
        rng = np.random.default_rng(0)
        r = circ_corr(a, b)
        null = np.array([circ_corr(a, rng.permutation(b)) for _ in range(n_perm)])
        rows.append(dict(analysis=lab, n=len(tl), r=r, p=(1 + np.sum(null >= r)) / (1 + n_perm), null_q95=float(np.quantile(null, .95)),
                         status="ok", tissues=";".join(tl)))
    P = pd.DataFrame(rows)
    P.to_csv(out / "a34_primary.csv", index=False)
    # exploratory: amplitude ranks
    mm = [(k, v) for k, v in MATCH.items() if k in tab.index and v in hum.index]
    amp = pd.DataFrame(dict(sheep_amplitude=[tab.rel_amplitude[k] for k, _ in mm], human_abs_g=[hum.abs_g[v] for _, v in mm]),
                       index=[k for k, _ in mm])
    amp = amp.dropna()
    rho = stats.spearmanr(amp.sheep_amplitude, amp.human_abs_g)
    amp.to_csv(out / "a34_amplitude.csv")
    pd.set_option("display.width", 250)
    print(P.drop(columns=["tissues"]).round(4).to_string(index=False))
    print(f"exploratory Spearman(sheep amplitude, human timed gain) = {rho.statistic:.3f} (p {rho.pvalue:.3g}, n {len(amp)})")


if __name__ == "__main__":
    main()
