"""A200 (registered in `paper/docs/registered_experiments.md`, commit b76dd7c): registered rerun of V25.5 on GSE233242 (42 breast
resections, tumour + adjacent normal, recorded hour 08-18 h). Is the per-patient internal time clock-specific?

Estimator unchanged (`results/figure_1/scripts/windows/compute_breast_timed.py`: two-compartment host, A25 estimate, flip by the
hospital prior, 30 M reads, pooled human waveform).
H200.1: Spearman(tangent, recorded hour) above the q95 of 200 expression-matched pseudo-clock sets run through the same pipeline
        (each clock gene replaced by one of its 20 nearest non-clock-like, non-G2M genes in mean and SD of log2(TPM + 1); the
        pseudo genes take the clock genes' names, so they inherit the waveform's amplitudes and phases).
H200.2: median |error| with the estimates' circular mean set to 13:00 (no recorded time used) below 2,000 shuffles of the hours.
Holm over two. Recorded hospital times of a living cohort are not controlled access; outputs are aggregates anyway.

    python score_a200.py [--npseudo 200] [--nperm 2000]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path[:0] = [str(PAPER), str(PAPER / "results/figure_1/scripts/windows")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
import compute_breast_timed as V  # noqa: E402
import _internal_time as IT  # noqa: E402
import importlib.util  # noqa: E402
from common.datasets.breast_timed_gse233242 import build as B  # noqa: E402

_s = importlib.util.spec_from_file_location("cancer_shared", PAPER / "results/figure_5/scripts/cancer/_shared.py")
_cancer = importlib.util.module_from_spec(_s); _s.loader.exec_module(_cancer)
hallmark_g2m = _cancer.hallmark_g2m
OUT = HERE.parents[1] / "output"
SEED, N_READS, PRIOR = 20261103, 30e6, 13.0
CLOCK = list(cs.gtex().genes)


def tables():
    """(TPM genes x libraries, library key, donors) as in the builder's primary join."""
    t = pd.read_csv(B.TRUTH)
    t = t[t.note.fillna("") == ""]
    key = pd.concat([pd.DataFrame(dict(patient=t.patient_id, compartment=c, SL=t[f"tpm_column_{c}"]))
                     for c in ("normal", "tumour")], ignore_index=True)
    x = pd.read_csv(B.TPM).drop_duplicates("gene_name").set_index("gene_name")
    _, donors = B.build()
    return x, key, donors


def pseudobulk(x, key, genes):
    """Builder counts for `genes`, columns renamed to the clock genes (genes[i] stands in for CLOCK[i])."""
    rl = x[key.SL].mul(x.transcriptLen, axis=0)
    cnt = (rl / rl.sum() * N_READS).loc[genes].T.round()
    cnt.columns = CLOCK
    cnt.index = range(len(key))
    return cs.from_bulk(cnt, key, donor="patient", compartment="compartment", depth=pd.Series(N_READS, index=cnt.index))


def pseudo_sets(x, key, n, rng):
    L = np.log2(x[key.SL] + 1)
    excl = set(hallmark_g2m())
    pool = [g for g in L.index if canonical(g) not in CLOCK_LIKE and g not in excl and g not in CLOCK and L.loc[g].mean() > 1]
    m, s = L.loc[pool].mean(1), L.loc[pool].std(1)
    P = np.column_stack([(m - m.mean()) / m.std(), (s - s.mean()) / s.std()])
    T = np.column_stack([(L.loc[CLOCK].mean(1) - m.mean()) / m.std(), (L.loc[CLOCK].std(1) - s.mean()) / s.std()])
    sets = []
    for _ in range(n):
        used, out = set(), {}
        for i in rng.permutation(len(CLOCK)):
            near = [j for j in np.argsort(((P - T[i]) ** 2).sum(1))[:40] if j not in used][:20]
            j = int(rng.choice(near)); used.add(j); out[i] = pool[j]
        sets.append([out[i] for i in range(len(CLOCK))])
    return sets


def readout(pb, donors, drop=None):
    """(Spearman of the tangent with the hour, median |error| with the 13:00 rule, estimates, truth)."""
    ds = cs.Dataset(pb, anchors=V.COMPS, name="a200", waveform="human", drop_genes=drop)
    e, host = V.two_host_estimate(ds)
    truth = donors.resection_h.reindex(host.index).to_numpy(float) % 24
    est = e["donor_day_h"].to_numpy(float)
    rho = spearmanr(e["tangent_day"], truth)[0]
    return rho, err13(est, truth), est, truth


def err13(est, truth):
    off = IT.cdiff(PRIOR, (np.angle(np.exp(1j * IT.W * est).mean()) / IT.W) % 24)
    return float(np.median(np.abs(IT.cdiff(est + off, truth))))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--npseudo", type=int, default=200); ap.add_argument("--nperm", type=int, default=2000)
    a = ap.parse_args()
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    x, key, donors = tables()
    rho, e13, est, truth = readout(pseudobulk(x, key, CLOCK), donors)
    rng = np.random.default_rng(SEED)
    null13 = np.array([err13(est, rng.permutation(truth)) for _ in range(a.nperm)])
    p2 = (1 + np.sum(null13 <= e13)) / (1 + a.nperm)
    ps, fails = [], 0
    for s in pseudo_sets(x, key, a.npseudo, np.random.default_rng([SEED, 1])):
        try:
            r, e, _, _ = readout(pseudobulk(x, key, s), donors)
            ps.append((r, e))
        except Exception:
            fails += 1
    ps = np.array(ps)
    p1 = (1 + np.sum(ps[:, 0] >= rho)) / (1 + len(ps))
    fail_share = fails / a.npseudo
    pa = np.array([p1, p2])
    o = np.argsort(pa)
    holm = np.empty(2); holm[o] = np.minimum(1, np.maximum.accumulate(np.array([2, 1]) * pa[o]))
    tests = pd.DataFrame([
        dict(test="H200.1 Spearman(tangent, hour) vs pseudo-clock sets", value=rho, pseudo_median=float(np.median(ps[:, 0])),
             pseudo_q95=float(np.quantile(ps[:, 0], .95)), pseudo_percentile=float(np.mean(ps[:, 0] < rho) * 100), p=p1,
             pseudo_fits=len(ps), pseudo_failed=fails, informative=fail_share <= 0.2),
        dict(test="H200.2 median |error|, 13:00 rule (h)", value=e13, null_median=float(np.median(null13)), p=p2,
             pseudo_median=float(np.median(ps[:, 1])), pseudo_percentile=float(np.mean(ps[:, 1] > e13) * 100))])
    tests["p_holm"] = holm
    tests["verdict"] = np.where((tests.test.str.startswith("H200.1")) & (fail_share > 0.2), "uninformative",
                                np.where(tests.p_holm < .05, "pass", "fail"))
    tests.to_csv(OUT / "tests.csv", index=False)

    pb = pseudobulk(x, key, CLOCK)
    sec = []
    for lab, drop in (("without PER1", ["PER1"]), ("without PER1 and NR1D1", ["PER1", "NR1D1"])):
        r, e, _, _ = readout(pb, donors, drop)
        sec.append(dict(analysis=lab, spearman=r, err13_h=e))
    row, *_ = V.run(pb, donors, n_perm=0)
    sec.append(dict(analysis="normal alone (tangent of host)", spearman=row["spearman_normal_alone"]))
    sec.append(dict(analysis="tumour alone (tangent of host)", spearman=row["spearman_tumour_alone"]))
    row_b, *_ = V.run(pb, donors, waveform="Breast - Mammary Tissue", n_perm=0)
    sec.append(dict(analysis="breast waveform", spearman=row_b["spearman_t"]))
    sec = pd.DataFrame(sec)
    sec.to_csv(OUT / "secondary.csv", index=False)
    pd.set_option("display.width", 250)
    print(tests.round(4).to_string(index=False)); print(sec.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
