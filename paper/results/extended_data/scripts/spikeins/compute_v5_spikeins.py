"""V5: calibration and power of the reference-free contrasts (offset difference, log amplitude ratio) on real biopsies without times.

Port of the archive `scripts/v5_spikein.py`, design unchanged. Study SRP445031 (Human Breast Cell Atlas, normal breast, 20 biopsies),
epithelial compartment. The compartment is split into two binomial halves A and B (exchangeable, true contrast 0), so the coverage of 0
calibrates the intervals at the study's real n, depth and noise. A's coupled clock is then damped (x0.7, 0.5, 0.3) or shifted (2 h, 4 h)
by thinning along the host state estimated without A and B; the expected change is diluted by the reliability r of that host state,
E[z_A'] = z_A (1 + (q - 1) r), q = scale exp(i shift). "decouple": A's resections shuffled. 40 replicates per manipulation, seeds as in the
archive (`default_rng(0)` draws one seed per task, effect-major). Contrasts: `Dataset.joint_contrasts` (200 bootstraps).
Readouts: detection (95 % interval excludes 0), coverage of the expected change, by amplitude and offset.

    python compute_v5_spikeins.py [--subsample N] [--reps 40] [--jobs 4] [--umi U]
`--subsample N`: check mode (N replicates per effect, 20 bootstraps).
Outputs: paper/results/extended_data/output/spikeins/v5_{runs,summary}.csv, compare_archive.csv
Archive: scCircadianCancer/outputs/v5_spikein/SRP445031_epithelial.csv (input pseudobulk: outputs/pilot/SRP445031_pseudobulk.csv)
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402
from common import cli  # noqa: E402
from common.paths import ARCHIVE, out_dir  # noqa: E402

FOLDER = "spikeins"
STUDY, COMP, TISSUE, REF = "SRP445031", "epithelial", "Breast", "fibroblast"
PB = ARCHIVE / "outputs" / "pilot" / f"{STUDY}_pseudobulk.csv"
ARCHIVE_RUN = ARCHIVE / "outputs" / "v5_spikein" / f"{STUDY}_{COMP}.csv"
EFFECTS = [("null", 1.0, 0.0), ("damp_0.7", 0.7, 0.0), ("damp_0.5", 0.5, 0.0), ("damp_0.3", 0.3, 0.0),
           ("shift_2h", 1.0, 2.0), ("shift_4h", 1.0, 4.0), ("decouple", None, None)]
WF = cs.waveform.gtex(TISSUE)
GENES = list(WF.genes)


def split(pb: pd.DataFrame, k: str, rng) -> pd.DataFrame:
    """Compartment k -> k:A and k:B, binomial halves of its counts and depth."""
    i = pb.compartment == k
    a, b = pb[i].copy(), pb[i].copy()
    ca = rng.binomial(pb.loc[i, GENES].to_numpy().astype(np.int64), 0.5)
    a[GENES], b[GENES] = ca, pb.loc[i, GENES].to_numpy() - ca
    a["depth"] = b["depth"] = pb.loc[i, "depth"] / 2
    a["compartment"], b["compartment"] = f"{k}:A", f"{k}:B"
    return pd.concat([pb[~i], a, b], ignore_index=True)


def dataset(t: pd.DataFrame, A: str, B: str) -> cs.Dataset:
    """Anchors = every compartment but A and B (A and B are targets read on that host); joint model referenced on fibroblast."""
    anchors = [REF] + sorted(set(t.compartment) - {A, B, REF})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = cs.Dataset(t, anchors=anchors, waveform=WF, name=STUDY)
    assert ds.study.ref == REF
    return ds


def one(pb, effect, scale, shift, rep, seed, n_boot):
    cs.config.limit_threads(1)
    rng = np.random.default_rng(seed)
    k = COMP
    t = split(pb, k, rng)
    A, B = f"{k}:A", f"{k}:B"
    st = dataset(t, A, B).study
    host = st.fit(exclude=[A, B])
    r = 1.0 - np.trace(host.f_cov.mean(0)) / np.trace(host.sigma_f)  # reliability of the host state
    i = np.flatnonzero(st.compartment == A)
    if effect == "decouple":
        t.loc[i, "donor"] = t.donor.to_numpy()[rng.permutation(i)]
    elif effect != "null":  # same reference as the full fit, hence the same frame
        ft = st.fit()
        z_a = ft.z[list(ft.compartments).index(A)]
        pos = {s: j for j, s in enumerate(host.resections)}
        f = host.f[[pos[s] for s in st.resection[i]]]
        w = (f[:, 0] + 1j * f[:, 1]) * z_a  # the host state seen through A's loading
        coupled = np.column_stack([w.real, w.imag])
        t.loc[i, GENES] = cs.simulate.thin_clock(t.loc[i, GENES].to_numpy(), WF.X,
                                                 cs.simulate.clock_change(coupled, scale, shift), rng)
    ds2 = dataset(t, A, B)
    c = ds2.joint_contrasts([(A, B)], reference=REF, n_boot=n_boot, seed=seed).iloc[0].to_dict()
    q = None if scale is None else scale * np.exp(1j * shift / HOURS)
    expect = np.nan if q is None else 1 + (q - 1) * r
    c.update(study=STUDY, compartment=k, effect=effect, rep=rep, reliability=r,
             expect_log_amp_ratio=np.log(abs(expect)) if q is not None else np.nan,
             expect_offset_diff_h=np.angle(expect) * HOURS if q is not None else np.nan,
             coupling_A=float(ds2.joint(REF).fit().table().coupling[A]))
    return c


def readouts(res: pd.DataFrame) -> pd.DataFrame:
    res = res.copy()
    res["detect_amp"] = (res.log_amp_ratio_hi < 0) | (res.log_amp_ratio_lo > 0)
    res["detect_off"] = (res.offset_diff_hi < 0) | (res.offset_diff_lo > 0)
    res["cover_amp"] = (res.log_amp_ratio_lo <= res.expect_log_amp_ratio) & (res.expect_log_amp_ratio <= res.log_amp_ratio_hi)
    res["cover_off"] = (res.offset_diff_lo <= res.expect_offset_diff_h) & (res.expect_offset_diff_h <= res.offset_diff_hi)
    s = res.groupby("effect", sort=False).agg(
        n=("rep", "size"), detect_amp=("detect_amp", "mean"), detect_off=("detect_off", "mean"), cover_amp=("cover_amp", "mean"),
        cover_off=("cover_off", "mean"), log_amp_ratio=("log_amp_ratio", "mean"), expect_log_amp_ratio=("expect_log_amp_ratio", "mean"),
        offset_diff_h=("offset_diff_h", "mean"), expect_offset_diff_h=("expect_offset_diff_h", "mean"),
        coupling_A=("coupling_A", "mean"), reliability=("reliability", "mean"))
    return s


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--reps", type=int, default=40)
    ap.add_argument("--umi", type=float, default=None)
    args = ap.parse_args()
    opt = cli.setup(args)
    quick = opt["quick"]
    reps = args.subsample if quick else args.reps
    n_boot = 20 if quick else 200
    out = out_dir(FOLDER, quick)
    pb = pd.read_csv(PB).drop(columns=["donor", "kind"]).rename(columns={"resection": "donor"})  # a biopsy is the unit (archive "resection")
    assert not args.umi, "--umi (depth thinning) not ported: archive umi1500/umi3000 runs are not part of the V5 headline"
    rng = np.random.default_rng(0)  # archive order: effect-major, one seed per task
    tasks = [(e, s, h, rep, int(rng.integers(2**31))) for e, s, h in EFFECTS for rep in range(40)]
    tasks = [t for t in tasks if t[3] < reps]
    res = pd.DataFrame(Parallel(n_jobs=args.jobs, verbose=5)(delayed(one)(pb, *t, n_boot) for t in tasks))
    res.to_csv(out / "v5_runs.csv", index=False)
    new = readouts(res)
    new.to_csv(out / "v5_summary.csv")
    print(new.round(3).to_string(), flush=True)
    old_runs = pd.read_csv(ARCHIVE_RUN, keep_default_na=False, na_values=["", "nan", "NaN"])  # effect "null" is not a missing value
    old = readouts(old_runs)
    rows = []
    for eff in new.index:
        for m in new.columns:
            if m == "n":
                continue
            o, n_ = old.loc[eff, m], new.loc[eff, m]
            rows.append(dict(effect=eff, metric=m, archive=o, new=n_, difference=n_ - o))
    cmp = pd.DataFrame(rows)
    cmp.to_csv(out / "compare_archive.csv", index=False)
    # paired replicates (same seeds): how far do the single-run contrasts move
    m = res.merge(old_runs, on=["effect", "rep"], suffixes=("_new", "_old"))
    for c in ("log_amp_ratio", "offset_diff_h", "log_amp_ratio_lo", "log_amp_ratio_hi", "offset_diff_lo", "offset_diff_hi", "reliability"):
        m[f"d_{c}"] = m[f"{c}_new"] - m[f"{c}_old"]
    p = m.groupby("effect", sort=False)[[c for c in m if c.startswith("d_")]].agg(lambda x: x.abs().max())
    p.to_csv(out / "compare_archive_paired_maxabs.csv")
    print(cmp.round(4).to_string(), "\n", p.round(6).to_string())


if __name__ == "__main__":
    main()
