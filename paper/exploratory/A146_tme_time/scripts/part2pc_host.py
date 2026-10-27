"""A146 part 2 (pan-cancer) blind preparation, step 2: stromal host, gate and A25 internal time per cohort (tumour samples).

Blind: uses only the 14 clock genes of the stromal compartments (fibroblast, endothelial, mural of the tumour sample);
no program score is computed or read. Stromal anchors = the learned A17 anchors of the cohort present among the tumour stromal
compartments (A138 `stromal_anchors`: all tumour stromal compartments if fewer than 2 learned). Patients with a host =
patients with >= 2 of those compartments. Gate: `Dataset.gate` (1,000 permutations, seed 146); passes = both p < 0.05
(`Gate.passes`). Internal time: A25 window of the host (`_internal_time.estimate`), flip resolved by the hospital prior
(centre nearer 13:00): `donor_day_h`. Three or more anchors: `_internal_time.host_estimate`; exactly two: frame shift read from
the complex cross-correlation of the two compartments (as figure_1 `compute_breast_timed.two_host_estimate`).

    python part2pc_host.py [--jobs 2]
Outputs: output/part2pc_cache/host_cohorts.csv (one row per cohort), output/part2pc_cache/host_times.csv (per-patient
internal time, hours; used by the permutation step, never reported per patient).
"""
from __future__ import annotations
import sys, time, argparse, importlib.util
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed

REPO = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(REPO / "paper/exploratory/A138_cancer_estimability/scripts"),
                str(REPO / "paper/results/figure_1/scripts/windows")]
import q3_lib as q  # noqa: E402
import _internal_time as IT  # noqa: E402
import clocksync as cs  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output" / "part2pc_cache"
SEED, N_PERM = 146, 1000


def circ_sd_h(h) -> float:
    R = abs(np.exp(1j * np.asarray(h) / IT.HOURS).mean())
    return float(np.sqrt(-2 * np.log(R)) * IT.HOURS) if R > 0 else np.nan


def two_anchor_estimate(ds: cs.Dataset, comps: list):
    ft = ds.fit(n_boot=0, gate=False)
    d = ds.data
    b = {k: pd.DataFrame(d.beta[d.compartment == k], index=d.resection[d.compartment == k]) for k in comps}
    common = b[comps[0]].index.intersection(b[comps[1]].index)
    o = np.angle(cs.stats.ccorr(b[comps[1]].loc[common].to_numpy(), b[comps[0]].loc[common].to_numpy()))
    return IT.estimate(ft.host[["x", "y"]].to_numpy(), ft.sigma_f, o / 2, None, ft.host.index)


def one(study: str) -> tuple[dict, pd.DataFrame]:
    from clocksync.config import limit_threads
    limit_threads(1)
    pb, mal = q.study_table(study)
    spec = importlib.util.spec_from_file_location("b", REPO / "paper/common/datasets/cancer_sc/build.py")
    row = dict(study=study, malignant_compartment=mal)
    tum = pb[pb.compartment.str.startswith("tumour:")]
    row["n_patients_tumour"] = int(tum.donor.nunique())
    anc, src = q.stromal_anchors(study, list(tum.compartment.unique()))
    row.update(anchors=";".join(a.split(":")[1] for a in anc), n_anchors=len(anc), anchor_source=src)
    out = pd.DataFrame()
    try:
        t = tum[tum.compartment.isin(anc)]
        n = t.groupby("donor").compartment.nunique()
        t = t[t.donor.isin(n[n >= 2].index)]
        row["n_patients_stromal_host"] = int(t.donor.nunique())
        ds = cs.Dataset(t, anchors=anc, name=study)
        g = ds.gate(n_perm=N_PERM, seed=SEED)
        row.update(gate_clock_share=g.clock_share, gate_share_null=g.share_null_mean, gate_p_clock_share=g.p_acrophase,
                   gate_p_pairing=g.p_pairing, gate_passes=bool(g.passes))
        e = IT.host_estimate(ds, anc)[0] if len(anc) >= 3 else two_anchor_estimate(ds, anc)
        h = e["donor_day_h"]
        row.update(n_patients_internal_time=int(h.notna().sum()), time_circ_sd_h=circ_sd_h(h.dropna()),
                   time_centre_h=float(e["c_day"]), window_w_h=float(e["w"]), window_r=float(e["r"]))
        out = pd.DataFrame(dict(study=study, patient=h.index, internal_h=h.to_numpy()))
    except Exception as ex:  # noqa: BLE001
        row["error"] = f"{type(ex).__name__}: {ex}"[:200]
    print(study, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    return row, out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--jobs", type=int, default=2); a = ap.parse_args()
    spec = importlib.util.spec_from_file_location("b", REPO / "paper/common/datasets/cancer_sc/build.py")
    sys.path[:0] = [str(REPO / "paper"), str(REPO / "paper/common/datasets")]
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    res = Parallel(n_jobs=a.jobs)(delayed(one)(s) for s in m.STUDIES)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([r for r, _ in res]).to_csv(OUT / "host_cohorts.csv", index=False)
    pd.concat([t for _, t in res], ignore_index=True).to_csv(OUT / "host_times.csv", index=False)
