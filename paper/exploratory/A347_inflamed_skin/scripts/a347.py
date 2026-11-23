"""A347 (registered in `paper/docs/registered_experiments.md`): is the lesional clock delay of inflamed skin specific to
psoriasis, or shared by atopic dermatitis (AD)?  17 gated discovery series (A53 gate), unit = patient, lesional minus
non-lesional, readouts: the 14-gene waveform (as A54) and the core loop (`a336_common.core`).

  --stage power   swap nulls only: SD of the pooled swap-null offset per disease x readout (10,000 swaps), MDE = 2.49 SD,
                  expected half-width of D = 1.645 sqrt(SD_PSO^2 + SD_AD^2), decidability flags. No observed offset.
  --stage score   (needs the user's go-ahead) pooled offsets, one-sided swap p, companion calibration q (A308 form),
                  D = PSO - AD with a patient bootstrap, per-series descriptive offsets.
  --smoke         score code path on data whose lesional labels are randomly swapped within patient (no real value).

Weights: 1 / se^2 of the series' variance ratio r (bootstrap SE, `paired.boot_se`). The power stage takes them from the
stored A53 table (14-gene r, same weights for both readouts: it computes nothing on the real labels); the score stage
recomputes them per readout on its own waveform.
Seeds: companions [347, disease, family] (disease 0 psoriasis, 1 AD; family 0 matched, 1 responsive, 2 broad), swaps
[347, 99, disease] (one stream, series in sorted order, the first 2,000 rows serve the companion calibration, all
10,000 the clock's swap p and the power stage).

    python a347.py --stage power
    python a347.py --stage score [--smoke]
Outputs (`../output/`): power/power.csv; score/*.csv (smoke: smoke/).
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A336_two_readouts/scripts"))
sys.path.insert(0, str(PAPER / "common/datasets/_bulk_common"))
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common.paths import DATASETS  # noqa: E402
import a336_common as C  # noqa: E402
import bulk_tables as T  # noqa: E402
import importlib.util  # noqa: E402

OUT = HERE.parents[1] / "output"
A53 = PAPER / "results/figure_5/output/psoriasis_details/a53_skin_per_series.csv"
DISEASES = ["psoriasis", "atopic_dermatitis"]
FAMILIES = ["matched", "responsive", "responsive_broad"]
N_SWAP, N_CAL, N_SETS, N_BOOT_D = 10000, 2000, 200, 2000
LEVELS = ("nonlesional", "lesional")
WF14 = cs.gtex()
READOUTS = {"14": WF14, "core": C.core(WF14)}


def builder():
    f = DATASETS / "skin_paired_discovery" / "build.py"
    spec = importlib.util.spec_from_file_location("build_skin", f)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def gated(disease):
    t = pd.read_csv(A53)
    t = t[(t.disease == disease) & t.gated.astype(bool)]
    return sorted(t.series), (1 / t.set_index("series").se ** 2).to_dict()


def align(x):
    case = x[x.condition == LEVELS[1]].drop_duplicates("unit").set_index("unit")
    ref = x[x.condition == LEVELS[0]].drop_duplicates("unit").set_index("unit")
    u = sorted(set(case.index) & set(ref.index))
    return case.loc[u], ref.loc[u]


def pairs_14(disease, b):
    """Series -> (case, ref, counts) of the gated series, 14-gene table (builder)."""
    s, _ = b.build()
    keys, _ = gated(disease)
    out = {}
    for k in keys:
        x = s[(s.cohort == f"{k}|{disease}")]
        case, ref = align(x)
        out[k] = (case, ref, bool((x.value_type == "counts").all()))
    return out


def pairs_all(disease, b, genes_extra=True):
    """Series -> (case, ref, counts) with every gene column (companions need them), clock genes in panel symbols."""
    keys, _ = gated(disease)
    out = {}
    for k in keys:
        d = T.read_all_genes(b.SRC / f"{k}.parquet", T.META)
        d = d[d.disease == disease]
        gc = T.clock_columns(list(d.columns))
        d = d.rename(columns=dict(zip(gc, T.GENES)))
        d = d.rename(columns={"patient": "unit", "status": "condition"}).assign(unit=lambda z: z["unit"].astype(str))
        case, ref = align(d)
        out[k] = (case.dropna(axis=1, how="all"), ref.dropna(axis=1, how="all"), bool((d.value_type == "counts").all()))
    return out


def smoke_swap(data, seed=0):
    """Randomly swap the lesional / non-lesional rows of every patient (no real label value is read afterwards)."""
    rng = np.random.default_rng(seed)
    out = {}
    for k, (case, ref, c) in data.items():
        f = rng.random(len(case)) < 0.5
        a, r = case.copy(), ref.copy()
        a.iloc[f], r.iloc[f] = ref.iloc[f].to_numpy(), case.iloc[f].to_numpy()
        out[k] = (a, r, c)
    return out


# ---- statistics ------------------------------------------------------------------------------------------------

def complex_pairs(case, ref, counts, wf):
    """Centred complex states (lesional a, non-lesional b) and measurement variances, one projector (A308)."""
    ba, va, bb, vb = P.cohort_states(case, ref, counts, wf)
    ba, bb = ba - ba.mean(0), bb - bb.mean(0)
    return ba[:, 0] + 1j * ba[:, 1], va, bb[:, 0] + 1j * bb[:, 1], vb


def ccorr_draws(a, va, b, vb, F):
    """`stats.ccorr` of (lesional, non-lesional) under each row of swaps F (draws x persons), vectorised."""
    A_, B_ = np.where(F, b, a), np.where(F, a, b)
    VA, VB = np.where(F, vb, va), np.where(F, va, vb)
    A_, B_ = A_ - A_.mean(1, keepdims=True), B_ - B_.mean(1, keepdims=True)
    k = a.size - 1
    ea = np.sum(abs(A_) ** 2, 1) / k - VA.mean(1)
    eb = np.sum(abs(B_) ** 2, 1) / k - VB.mean(1)
    return np.sum(A_ * np.conj(B_), 1) / k / np.sqrt(np.maximum(ea, 1e-12) * np.maximum(eb, 1e-12))


def offsets(pairs, w, F, chunk=2000):
    """Pooled offset [h] per draw (rows of F[k]): arg sum_k w_k ccorr_k."""
    n = len(next(iter(F.values())))
    out = np.empty(n)
    for i in range(0, n, chunk):
        z = sum(w[k] * ccorr_draws(*pairs[k], F[k][i:i + chunk]) for k in pairs)
        out[i:i + chunk] = np.angle(z) * P.HOURS
    return out


def wrap(h):
    return (np.asarray(h, float) + 12.0) % 24.0 - 12.0


def sd0(h):
    """Circular-safe SD of a null centred on 0: root mean square of the offsets wrapped to [-12, 12) hours."""
    return float(np.sqrt(np.mean(wrap(h) ** 2)))


def wsd(h):
    """Placebo SD (A308): RMS of the offsets wrapped around their circular mean."""
    mu = np.angle(np.mean(np.exp(1j * np.asarray(h) * np.pi / 12))) * 12 / np.pi
    return float(np.sqrt(np.mean(wrap(h - mu) ** 2)))


def swap_matrix(disease_idx, sizes, n):
    """Series -> (n, persons) bool swap draws from the [347, 99, disease] stream (series in the order of `sizes`)."""
    rng = np.random.default_rng([347, 99, disease_idx])
    return {k: rng.random((n, m)) < 0.5 for k, m in sizes.items()}


# ---- power -----------------------------------------------------------------------------------------------------

def power():
    t0 = time.time()
    b = builder()
    rows, units = [], []
    for d, dis in enumerate(DISEASES):
        data = pairs_14(dis, b)
        _, w = gated(dis)
        sizes = {k: len(v[0]) for k, v in data.items()}
        units.append(dict(disease=dis, n_series=len(sizes), n_pairs=sum(sizes.values()), **{f"n_{k}": v for k, v in sizes.items()}))
        F = swap_matrix(d, sizes, N_SWAP)
        for name, wf in READOUTS.items():
            pr = {k: complex_pairs(*data[k], wf) for k in data}
            null = offsets(pr, w, F)
            rows.append(dict(disease=dis, readout=name, n_series=len(sizes), n_pairs=sum(sizes.values()), n_swaps=N_SWAP,
                             null_mean_h=float(wrap(null).mean()), sd_h=sd0(null), circ_sd_h=float(
                                 np.sqrt(-2 * np.log(abs(np.mean(np.exp(1j * null * np.pi / 12))))) * 12 / np.pi)))
        print(dis, "units", units[-1], flush=True)
    T_ = pd.DataFrame(rows)
    T_["mde_h"] = 2.49 * T_.sd_h
    for name in READOUTS:
        sd = T_[T_.readout == name].set_index("disease").sd_h
        hw = 1.645 * np.sqrt(sd.psoriasis ** 2 + sd.atopic_dermatitis ** 2)
        T_.loc[T_.readout == name, "D_halfwidth90_h"] = hw
    T_["ad_mde_le_0.75h"] = (T_.disease == "atopic_dermatitis") & (T_.mde_h <= 0.75)
    T_["D_halfwidth_lt_0.75h"] = T_.D_halfwidth90_h < 0.75
    (OUT / "power").mkdir(parents=True, exist_ok=True)
    T_.to_csv(OUT / "power/power.csv", index=False)
    pd.DataFrame(units).to_csv(OUT / "power/unit_counts.csv", index=False)
    print(T_.round(3).to_string(index=False))
    print(f"runtime {time.time() - t0:.0f} s")


# ---- score -----------------------------------------------------------------------------------------------------

def weights_of(data, wf, d, readout, n_boot):
    """1 / se^2 of the variance ratio r of each series on the readout waveform (bootstrap SE, A54 estimator)."""
    w = {}
    for i, (k, (case, ref, c)) in enumerate(data.items()):
        brng = np.random.default_rng([347, d, 9, i, int(readout == "core")])
        w[k] = 1 / P.boot_se(*P.cohort_states(case, ref, c, wf), n_boot=n_boot, rng=brng) ** 2
    return w


def boot_offset(pairs, w, rng, n):
    """Patient bootstrap within series: pooled offset [h] per replicate (weights fixed)."""
    out = np.empty(n)
    for r in range(n):
        z = 0
        for k, (a, va, b, vb) in pairs.items():
            i = rng.integers(0, len(a), len(a))
            z += w[k] * ccorr_draws(a[i], va[i], b[i], vb[i], np.zeros((1, len(a)), bool))[0]
        out[r] = np.angle(z) * P.HOURS
    return out


def score(smoke=False):
    t0 = time.time()
    n_swap, n_cal, n_sets, n_boot, n_bw = (500, 200, 20, 100, 50) if smoke else (N_SWAP, N_CAL, N_SETS, N_BOOT_D, 1000)
    out = OUT / ("smoke" if smoke else "score")
    out.mkdir(parents=True, exist_ok=True)
    b = builder()
    means = b.gene_means()
    res, fam_rows, ser_rows, comp = {}, [], [], []
    for d, dis in enumerate(DISEASES):
        data = pairs_all(dis, b)
        if smoke:
            data = smoke_swap(data, seed=d)
        sizes = {k: len(v[0]) for k, v in data.items()}
        Ffull = swap_matrix(d, sizes, n_swap)
        Fcal = {k: np.vstack([np.zeros((1, sizes[k]), bool), Ffull[k][:n_cal]]) for k in data}
        Fid = {k: np.zeros((1, sizes[k]), bool) for k in data}
        for ro, wf in READOUTS.items():
            w = weights_of(data, wf, d, ro, n_bw)
            pr = {k: complex_pairs(*data[k], wf) for k in data}
            o_obs = offsets(pr, w, Fid)[0]
            null = offsets(pr, w, Ffull)
            p_delay = S.perm_p(null, o_obs, "less", tie=0.0)
            clock = offsets(pr, w, Fcal)
            sc = wsd(clock[1:])
            for k in data:
                ser_rows.append(dict(disease=dis, readout=ro, series=k, n_pairs=sizes[k], weight=w[k],
                                     offset_h=float(np.angle(ccorr_one(pr[k])) * P.HOURS)))
            for f, fam in enumerate(FAMILIES):
                crng = np.random.default_rng([347, d, f])
                per = {k: P.pseudo_waveforms(means[f"{k}|{dis}"], n_sets, crng, WF14, fam) for k in data}
                M = np.empty((n_sets, n_cal + 1))
                for j in range(n_sets):
                    wfj = {k: (C.restrict(per[k][j], base=WF14) if ro == "core" else per[k][j]) for k in data}
                    M[j] = offsets({k: complex_pairs(*data[k], wfj[k]) for k in data}, w, Fcal)
                sj = np.array([wsd(m[1:]) for m in M])

                def q_draws(Mx, idx=slice(None)):
                    m = Mx[idx].mean(0)
                    tau2 = np.maximum(0.0, Mx[idx].var(0, ddof=1) - np.mean(sj[idx] ** 2))
                    return norm.cdf((clock - m) / np.sqrt(tau2 + sc ** 2)), tau2

                qd, tau2 = q_draws(M)
                brng = np.random.default_rng([347, d, f, 7])
                qb = [q_draws(M, brng.integers(0, n_sets, n_sets))[0][0] for _ in range(2000 if not smoke else 100)]
                fam_rows.append(dict(disease=dis, readout=ro, family=fam, clock_offset_h=o_obs, p_delay=p_delay,
                                     n_swaps=n_swap, clock_placebo_sd_h=sc, sets_mean_h=float(M[:, 0].mean()),
                                     tau2_obs=float(tau2[0]), q=float(qd[0]), q_lo90=float(np.quantile(qb, .05)),
                                     q_hi90=float(np.quantile(qb, .95)), v347_pass_rate=float(np.mean(qd[1:] < 0.05)),
                                     v347_calibrated=bool(np.mean(qd[1:] < 0.05) <= 0.0646)))
            res[(dis, ro)] = (pr, w, o_obs)
            print(dis, ro, "done", f"{time.time() - t0:.0f} s", flush=True)
    rows = []
    for ro in READOUTS:
        rng = np.random.default_rng([347, 5, int(ro == "core")])
        bo = {dis: boot_offset(res[(dis, ro)][0], res[(dis, ro)][1], rng, n_boot) for dis in DISEASES}
        o = {dis: res[(dis, ro)][2] for dis in DISEASES}
        D = wrap(o["psoriasis"] - o["atopic_dermatitis"])
        Db = D + wrap(wrap(bo["psoriasis"] - bo["atopic_dermatitis"]) - D)
        lo, hi = np.quantile(Db, [0.05, 0.95])
        rows.append(dict(readout=ro, offset_pso_h=o["psoriasis"], offset_ad_h=o["atopic_dermatitis"], D_h=D,
                         D_lo90=lo, D_hi90=hi, n_boot=n_boot, H347_3=bool(hi < 0), H347_4=bool(lo > -0.75 and hi < 0.75)))
    pd.DataFrame(rows).to_csv(out / "difference.csv", index=False)
    pd.DataFrame(ser_rows).to_csv(out / "per_series.csv", index=False)
    fam = pd.DataFrame(fam_rows)
    fam.to_csv(out / "families.csv", index=False)
    hyp = []
    for ro in READOUTS:
        g = fam[fam.readout == ro].set_index(["disease", "family"])
        ad, ps = g.loc["atopic_dermatitis"], g.loc["psoriasis"]
        hyp.append(dict(readout=ro, H347_1=bool(ad.p_delay.iloc[0] < 0.05),
                        H347_2=bool(ad.loc["matched", "q"] < 0.05 and ad.loc["responsive", "q"] < 0.05),
                        P347=bool(ps.p_delay.iloc[0] < 0.05 and ps.loc["matched", "q"] < 0.05)))
    pd.DataFrame(hyp).to_csv(out / "hypotheses.csv", index=False)
    print(fam.round(3).to_string(index=False), "\n", pd.DataFrame(rows).round(3).to_string(index=False),
          "\n", pd.DataFrame(hyp).to_string(index=False), f"\nruntime {time.time() - t0:.0f} s")


def ccorr_one(p):
    a, va, b, vb = p
    return ccorr_draws(a, va, b, vb, np.zeros((1, len(a)), bool))[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--smoke", action="store_true", help="score on randomly swapped labels, reduced sizes")
    ap.add_argument("--go", action="store_true", help="real-label scoring: only with the user's go-ahead")
    a = ap.parse_args()
    if a.stage == "score" and not (a.smoke or a.go):
        sys.exit("score on the real labels needs --go (user's go-ahead); use --smoke to test the code path")
    limit_threads(4)
    if a.stage == "power":
        power()
    else:
        score(a.smoke)


if __name__ == "__main__":
    main()
