"""A123: the GTEx organ map from working-hours deaths only (registered in paper/docs/registered_experiments.md, WP4 note).

Donors: GTEx v10 fast deaths (Hardy 1-2) who died 08-18 h (`window`); tissues with >= 20 such donors (whole blood included).
Untimed estimates from the window donors only (times hidden): `clocksync` reference-free tissue offsets and relative amplitudes
(`Dataset.fit`, all tissues anchors, donor bootstrap; waveform fitted on the fast deaths who died OUTSIDE the window, as the
benchmark's cross-fit) and the per-donor internal time (A25 `host_estimate`, tangent coordinate with the 12 h flip resolved by
the hospital prior).
Truth (i): per-tissue cosinor of the 14 clock genes on the recorded time of ALL fast deaths (all hours), rank-1 complex SVD
(benchmark `csvd_estimates`). Truth (ii): the same harmonic regression on CHIRAL's published per-sample CT (naef-lab/CHIRAL commit
1f78a9a, `sample_metadata.txt`, all Hardy classes with a CT; whole blood has none), re-derived (the per-tissue tables are not local).
H123.1 circular correlation of centred offsets vs truth (i) > 0.7, tissue-permutation p < 0.05; H123.2 Pearson of centred log2
relative amplitudes > 0.7; H123.3 Spearman of the tangent coordinate with recorded time among window donors > 0.4, direction
correct (positive). Secondary: the same against truth (ii); CHIRAL (pychiral) run per tissue on the same window donors, per-donor
phases combined by the benchmark's consensus, tissue estimands by cosinor + SVD on the consensus times (comparator that assigns phases).
Controlled access: only per-tissue aggregates and summary statistics are written; per-donor times, phases and the tangent
coordinate stay in memory.

    python compute_working_hours.py [--subsample N] [--boot B] [--perm P]
`--subsample N`: check mode, the N tissues with most window donors, 10 bootstrap replicates, 200 permutations.
Outputs (exploratory/A123_working_hours/output/): organ_map.csv (per tissue), tests.csv, truth_i.csv, truth_ii.csv.
"""

from __future__ import annotations

import importlib.util
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_1/scripts/method"))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats as sps  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.benchmark import evaluate as ev  # noqa: E402
from common.benchmark.chiral import CHIRAL_GENES, run_chiral  # noqa: E402
from common.benchmark.data import in_window, log2_normalised  # noqa: E402
from common.gtex import FAST, LABELLED  # noqa: E402
from common.paths import DATA, exploratory_dir  # noqa: E402

_spec = importlib.util.spec_from_file_location("gtex_benchmark", PAPER / "results/figure_1/scripts/method/compute_gtex_benchmark.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)
sys.path.insert(0, str(PAPER / "results/figure_1/scripts/windows"))
import _internal_time as T  # noqa: E402

AMENDMENT = "A123_working_hours"
CLOCK, REF, W = B.CLOCK, B.REF, ev.W
META = DATA / "chiral_repo" / "raw" / "sample_metadata.txt"
MIN_WINDOW = 20
SEED = 123


def window_tissues(o: pd.DataFrame) -> pd.Series:
    """Window donors per tissue (fast deaths with time of death), tissues with >= MIN_WINDOW."""
    o = o[~o.tissue.str.startswith("Cells -")]  # cultured cells are not tissues
    f = o[o.DTHHRDY.isin(list(FAST)) & o.tod.notna()].drop_duplicates(["subject", "tissue"])
    f = f[in_window(f.tod.astype(float) % 24)]
    n = f.groupby("tissue").subject.nunique()
    return n[n >= MIN_WINDOW].sort_values(ascending=False)




def truth_ii_coef(a, o, ct: pd.Series, tissue: str) -> pd.Series | None:
    """Complex cosinor coefficient per clock gene on CHIRAL's CT, all samples of the tissue with a CT (one per donor)."""
    m = (o.tissue == tissue).to_numpy() & o.index.isin(ct.index)
    idx = np.flatnonzero(m)
    if len(idx) < 20:
        return None
    X = a.X[idx]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    df = pd.DataFrame(X, index=o.index[idx], columns=a.var_names).T.groupby(level=0).sum().T
    return ev.gene_rhythms(log2_normalised(df), ct.loc[df.index], CLOCK)


def centre(est: pd.DataFrame, k) -> tuple[pd.DataFrame, float, float]:
    """Centred estimate over tissues k, and the circular centre (h) and mean log2 amplitude removed."""
    e = est.loc[list(k)]
    la = np.log2(e.amplitude)
    c = np.angle(np.exp(1j * W * e.offset_h).mean()) / W
    return ev.centred(est, k), c, la.mean()


def compare(est: pd.DataFrame, truth: pd.DataFrame, n_perm: int, seed: int) -> dict:
    """H123.1 / H123.2 statistics of centred offsets and log2 amplitudes (tissues common to both)."""
    k = est.dropna(subset=["amplitude", "offset_h"]).index.intersection(truth.index)
    te = ev.tissue_errors(est, truth)
    rng = np.random.default_rng(seed)
    r_off = ev.circ_corr(te.offset_h, te.offset_h_true)
    r_amp = float(np.corrcoef(te.log2_amp, te.log2_amp_true)[0, 1])
    po, pa = [], []
    for _ in range(n_perm):
        p = rng.permutation(len(k))
        po.append(ev.circ_corr(te.offset_h, te.offset_h_true.to_numpy()[p]))
        pa.append(np.corrcoef(te.log2_amp, te.log2_amp_true.to_numpy()[p])[0, 1])
    boot_o, boot_a = [], []
    for _ in range(1000):
        i = rng.integers(0, len(k), len(k))
        boot_o.append(ev.circ_corr(te.offset_h.to_numpy()[i], te.offset_h_true.to_numpy()[i]))
        boot_a.append(np.corrcoef(te.log2_amp.to_numpy()[i], te.log2_amp_true.to_numpy()[i])[0, 1])
    lin = np.corrcoef(2 ** te.log2_amp, 2 ** te.log2_amp_true)[0, 1]
    return dict(n_tissues=len(k), circ_corr_offset=r_off, p_offset=(1 + np.sum(np.array(po) >= r_off)) / (1 + n_perm),
                circ_corr_offset_lo=np.nanquantile(boot_o, .025), circ_corr_offset_hi=np.nanquantile(boot_o, .975),
                median_abs_offset_err_h=float(te.offset_err_h.abs().median()),
                pearson_log2_amp=r_amp, p_amp=(1 + np.sum(np.array(pa) >= r_amp)) / (1 + n_perm),
                pearson_log2_amp_lo=np.nanquantile(boot_a, .025), pearson_log2_amp_hi=np.nanquantile(boot_a, .975),
                pearson_linear_amp=float(lin), median_abs_log2_amp_err=float(te.log2_amp_err.abs().median()))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--boot", type=int, default=None)
    ap.add_argument("--perm", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    n_boot = args.boot if args.boot is not None else 10 if quick else 200
    n_perm = args.perm if args.perm is not None else 200 if quick else 10000
    out = exploratory_dir(AMENDMENT, quick)
    t0 = time.time()

    a = ad.read_h5ad(LABELLED, backed="r")
    o = a.obs[["tissue", "subject", "DTHHRDY", "tod"]].copy()
    o["tissue"] = o.tissue.astype(str)
    cnt = window_tissues(o)
    tissues = list(cnt.index[:args.subsample]) if quick else list(cnt.index)
    if REF not in tissues:
        tissues = [REF] + tissues[:-1] if quick else tissues
    print(f"{len(tissues)} tissues with >= {MIN_WINDOW} window donors", flush=True)

    meta = pd.read_csv(META, sep=" ")
    meta["sample"] = meta.fullID.str.replace(".", "-", regex=False)
    ct = meta.set_index("sample").CT.dropna()

    Ci, Cii, rows, chiral_per, L_win, tod_of = {}, {}, [], {}, {}, {}
    for t in tissues:
        counts, tod = B.S.load([t])[t]
        tod = tod.dropna()
        counts = counts.loc[tod.index]
        win = tod.index[in_window(tod)]
        tod_of[t] = tod
        rows.append(B.clock_table({t: (counts, tod)}))
        if len(counts) >= 5:
            Ci[t] = ev.gene_rhythms(log2_normalised(counts), tod, CLOCK)
        Lw = log2_normalised(counts.loc[win])
        L_win[t] = Lw[[g for g in CLOCK if g in Lw.columns]]
        phase, failed = run_chiral(Lw, seed=0)
        chiral_per[t] = phase.dropna()
        c2 = truth_ii_coef(a, o, ct, t)
        if c2 is not None:
            Cii[t] = c2
        print(t, f"fast {len(counts)} window {len(win)} chiral failed starts {failed}", f"{time.time() - t0:.0f} s", flush=True)
    a.file.close()
    table = pd.concat(rows, ignore_index=True)
    window_donors = sorted({d for t in tissues for d in tod_of[t].index[in_window(tod_of[t])]})
    tod_all = pd.concat(list(tod_of.values())).groupby(level=0).first()

    truth_i = ev.csvd_tissues(pd.DataFrame(Ci), REF)
    truth_ii = ev.csvd_tissues(pd.DataFrame(Cii), REF) if REF in Cii else None
    truth_i.to_csv(out / "truth_i.csv")
    if truth_ii is not None:
        truth_ii.to_csv(out / "truth_ii.csv")

    # ---- untimed clocksync estimates from the window donors -----------------------------------------------------------
    wtab = table[table.donor.isin(window_donors) & in_window(table.tod)]
    wf = B.cross_fit_waveform(table, window_donors)  # fast deaths outside the window only
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est = B.clocksync_estimates(wtab, wf, n_boot, seed=SEED)
        est_nb = B.clocksync_estimates(wtab, cs.gtex(), 0)  # sensitivity: the all-donor GTEx waveform

        # ---- H123.3: per-donor internal time, tangent coordinate with the hospital prior -------------------------------
        anchors = sorted(wtab.compartment.unique())
        ds = cs.Dataset(wtab.drop(columns="tod"), anchors=anchors, waveform=wf)
        e, _ = T.host_estimate(ds, anchors)
    tan, tw = e["tangent_day"].dropna(), tod_all
    tw = tw.reindex(tan.index)
    rho, p3 = sps.spearmanr(tan.to_numpy(), tw.to_numpy())
    rng = np.random.default_rng(SEED)
    bs = []
    for _ in range(1000):
        i = rng.integers(0, len(tan), len(tan))
        bs.append(sps.spearmanr(tan.to_numpy()[i], tw.to_numpy()[i])[0])
    h3 = dict(n_donors=len(tan), spearman=rho, spearman_lo=np.nanquantile(bs, .025), spearman_hi=np.nanquantile(bs, .975),
              p_two_sided=p3, c_axis12=e["c_axis12"], c_day_h=e["c_day"], r=e["r"], w_h=e["w"])
    tr = T.true_window(tw.to_numpy())
    h3.update(true_window_w_h=tr["w_true"], true_window_r=tr["r_true"])
    # the same with the curvature-resolved flip (no prior), descriptive
    t2 = e["tangent"].dropna()
    h3["spearman_curvature_flip"] = float(sps.spearmanr(t2.to_numpy(), tod_all.reindex(t2.index).to_numpy())[0])
    # ---- CHIRAL comparator -------------------------------------------------------------------------------------------
    cons = B.consensus_times(chiral_per, REF) if REF in chiral_per and len(chiral_per[REF]) else None
    if cons is not None:
        Cc = {}
        for t in tissues:
            ix = L_win[t].index.intersection(cons.index)
            if len(ix) >= 5:
                Cc[t] = ev.gene_rhythms(L_win[t].loc[ix], cons.loc[ix], CLOCK)
        chiral = ev.csvd_tissues(pd.DataFrame(Cc), REF)
        trw = tod_all.reindex(cons.index)
        _, mae, _ = ev.align(cons, trw)
        ch_don = dict(consensus_n=len(cons), consensus_circ_corr=abs(ev.circ_corr(cons, trw)), consensus_mae_h=mae)
        pt = [abs(ev.circ_corr(p, tod_of[t].reindex(p.index))) for t, p in chiral_per.items() if len(p) >= 5]
        ch_don["median_tissue_circ_corr"] = float(np.median(pt))
    else:
        chiral, ch_don = None, {}

    # ---- scores --------------------------------------------------------------------------------------------------------
    tests = []
    for name, es in (("clocksync", est), ("clocksync_gtex_waveform", est_nb), ("chiral", chiral)):
        if es is None:
            continue
        for tn, tr_ in (("truth_i", truth_i), ("truth_ii", truth_ii)):
            if tr_ is None:
                continue
            tests.append(dict(method=name, truth=tn, **compare(es, tr_, n_perm, SEED)))
    tests = pd.DataFrame(tests)
    tests["H123.1"] = np.where((tests.circ_corr_offset > 0.7) & (tests.p_offset < 0.05), "pass", "fail")
    tests["H123.2"] = np.where(tests.pearson_log2_amp > 0.7, "pass", "fail")
    tests.to_csv(out / "tests.csv", index=False)
    h3["H123.3"] = "pass" if (h3["spearman"] > 0.4) else "fail"
    h3.update({f"chiral_{k}": v for k, v in ch_don.items()})
    pd.Series(h3).to_csv(out / "donor_time.csv", header=False)

    # ---- per-tissue deliverable (reference-free, centred over the tissues common to estimate and truth (i)) ---------------
    k = est.index.intersection(truth_i.index)
    ec, off_c, la_c = centre(est, k)
    d = pd.DataFrame(index=pd.Index(k, name="tissue"))
    d["n_window_donors"] = [int(cnt.get(t, np.nan)) for t in k]
    d["est_log2_amp"] = ec.log2_amp
    d["est_log2_amp_lo"] = est.log2_amp_lo.loc[k] - la_c
    d["est_log2_amp_hi"] = est.log2_amp_hi.loc[k] - la_c
    d["est_offset_h"] = ec.offset_h
    d["est_offset_lo"] = ev.wrap(est.offset_lo.loc[k] - off_c)
    d["est_offset_hi"] = ev.wrap(est.offset_hi.loc[k] - off_c)
    d["coupling"] = est.coupling.loc[k]
    ti = ev.centred(truth_i, k)
    d["truth_i_log2_amp"], d["truth_i_offset_h"] = ti.log2_amp, ti.offset_h
    if truth_ii is not None:
        k2 = k.intersection(truth_ii.index)
        t2c = ev.centred(truth_ii, k2)
        d["truth_ii_log2_amp"], d["truth_ii_offset_h"] = t2c.log2_amp, t2c.offset_h
    if chiral is not None:
        kc = k.intersection(chiral.index)
        cc = ev.centred(chiral, kc)
        d["chiral_log2_amp"], d["chiral_offset_h"] = cc.log2_amp, cc.offset_h
    d.round(4).to_csv(out / "organ_map.csv")

    pd.set_option("display.width", 250)
    print(tests.round(3).T.to_string())
    print(pd.Series(h3).to_string())


if __name__ == "__main__":
    main()
