"""A387: tumour compartments against the patient's normal tissue, complex-pooled pairing and full model.

Registered in paper/docs/registered_experiments.md (A387). Cohorts, tables and units of A346 (A227 tables; A149 units for
the malignant-normal epithelium pairing). Readouts: 14 genes (the A227 / A336 waveform) and the core loop.

  Part 1  complex-pooled pairing R_C = |sum_k sum_i z_T conj(z_N)| / sqrt(sum |z_T|^2 sum |z_N|^2) over cohorts k and
          patients i, states as A346 (`A346_tumour_stroma_core/output/cache/units.pkl`: rows 0 / 601 clock, then the
          200 matched, 200 narrow responsive and 200 broad responsive sets of each readout). Phase = arg of the pooled
          cross-moment, tumour minus normal (clocksync sign: > 0 = tumour advanced).
  Part 2  one `clocksync.Dataset` over all A227 cohorts (donor and stratum prefixed by cohort), anchors = the normal
          sample's fibroblasts, endothelium and mural cells; targets the tumour fibroblasts, endothelium and malignant
          cells, and the normal epithelium.

  power   (no clock target statistic computed)
          Part 1: patients and cohorts per compartment; within-cohort permutation null of R_C (shuffled pairings only,
          10,000, seeds [387, kind, 1]), MDE = 2.49 SD; the companions on the real pairing (clock rows excluded); 2,000
          placebo draws (normal rows permuted within chemistry stratum within cohort, seeds [387, kind, d]; cached for
          the score stage); analytic 90 % half-width of the pooled phase at a hypothetical R_C.
          Part 2: patients per target with a host; gate on the anchors only; power of the coupling at the real design
          by the blind simulator (`clocksync.simulate.like(..., loading=z)`: the target is never regressed on the
          host), c in a grid, |z| 1 and 0.5 (host gauge), null threshold from c = 0.02; 90 % half-width of the target's
          offset at c 0.5, |z| 1. Decidability as registered.            -> output/power/*.csv, output/cache/
  score   (not written before the power note is committed)

    python a387.py --stage power [--jobs 4] [--smoke]     (long: paper/common/run.sh, RUN_CAP_GB 30)
"""
from __future__ import annotations

import argparse
import pickle
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(PAPER))
for p in ("A227_tme_compartment_pairing/scripts", "A149_paired_tumour_normal/scripts", "A151_fibroblast_pairing/scripts",
          "A226_stroma_pairing_confirm/scripts", "A336_two_readouts/scripts"):
    sys.path.insert(0, str(EXP / p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a336_common as C  # noqa: E402
import compute_a227 as M  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import HOURS  # noqa: E402

OUT = HERE.parents[1] / "output"
A346_CACHE = EXP / "A346_tumour_stroma_core/output/cache/units.pkl"
BASE = C.WF
READOUTS = {"14": BASE, "core": C.core(BASE)}
N_SET, N_PERM, N_DRAW, N_SIM = 200, 10_000, 2000, 200
NB = 1 + 3 * N_SET
BLOCK = {"14": 0, "core": NB}
KCODE = {"fibroblast": 1, "endothelial": 2, "epithelium": 3}
PRIMARY = ("fibroblast", "endothelial")
MDE_K, MAX_MDE, BOUND_H = 2.49, 0.25, 1.5
ANCHORS = ["adjacent_normal:fibroblast", "adjacent_normal:endothelial", "adjacent_normal:mural"]
TARGETS = {"fibroblast": "tumour:fibroblast", "endothelial": "tumour:endothelial", "epithelium": "tumour:malignant"}
NORMAL = {"fibroblast": "adjacent_normal:fibroblast", "endothelial": "adjacent_normal:endothelial",
          "epithelium": "adjacent_normal:epithelial"}
C_GRID = (0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


# ---- Part 1: complex-pooled pairing -----------------------------------------------------------------------------

def _sums(za, zb):
    """Per row: cross-moment and the two energies of centred states (NaN rows contribute nothing)."""
    za, zb = za - za.mean(-1, keepdims=True), zb - zb.mean(-1, keepdims=True)
    ok = np.isfinite(za).all(-1) & np.isfinite(zb).all(-1)
    za, zb = np.where(ok[..., None], za, 0), np.where(ok[..., None], zb, 0)
    return (za * np.conj(zb)).sum(-1), (np.abs(za) ** 2).sum(-1), (np.abs(zb) ** 2).sum(-1)


def cpool(units, rows, rng=None):
    """R_C and pooled phase of the given rows; rng: one placebo draw (normal rows permuted within stratum)."""
    num, ea, eb = 0j, 0.0, 0.0
    for u in units:
        zb = u["zb"][rows]
        if rng is not None:
            idx = np.arange(u["n"])
            for st in np.unique(u["stratum"]):
                j = np.flatnonzero(u["stratum"] == st)
                idx[j] = j[rng.permutation(len(j))]
            zb = zb[:, idx]
        c, a, b = _sums(u["za"][rows], zb)
        num, ea, eb = num + c, ea + a, eb + b
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.abs(num) / np.sqrt(ea * eb), np.angle(num)


def perm_null(units, row, kind):
    """R_C of the clock row under within-cohort permutations of the normal rows (shuffled pairings only)."""
    num, ea, eb = np.zeros(N_PERM, complex), 0.0, 0.0
    for u in units:
        rng = np.random.default_rng([387, KCODE[kind], 1, sum(map(ord, u["study"]))])
        za, zb = u["za"][row] - u["za"][row].mean(), u["zb"][row] - u["zb"][row].mean()
        perms = np.argsort(rng.random((N_PERM, u["n"])), axis=1)
        num += (za[None, :] * np.conj(zb[perms])).sum(1)
        ea, eb = ea + np.sum(np.abs(za) ** 2), eb + np.sum(np.abs(zb) ** 2)
    return np.abs(num) / np.sqrt(ea * eb)


def phase_halfwidth_h(rho, n):
    """Analytic 90 % half-width (h) of the phase of a complex correlation rho over n pairs (large-sample,
    circular-normal states): SD = sqrt((1 - rho^2) / (2 n rho^2)) rad."""
    return 1.645 * np.sqrt((1 - rho ** 2) / (2 * n * rho ** 2)) * HOURS


def part1_power(smoke):
    with open(A346_CACHE, "rb") as f:
        U = pickle.load(f)
    cache = OUT / ("cache_smoke" if smoke else "cache")
    cache.mkdir(parents=True, exist_ok=True)
    rows = []
    for kind, units in U.items():
        if smoke:
            units = units[:3]
        n = int(sum(u["n"] for u in units))
        D = np.array([cpool(units, np.arange(2 * NB), np.random.default_rng([387, KCODE[kind], d]))[0]
                      for d in range(1, N_DRAW + 1)])
        comp = {rd: cpool(units, np.arange(BLOCK[rd] + 1, BLOCK[rd] + NB))[0] for rd in BLOCK}
        np.savez_compressed(cache / f"part1_draws_{kind}.npz", draws=D, comp14=comp["14"], compcore=comp["core"])
        for rd in BLOCK:
            nul = perm_null(units, BLOCK[rd], kind)
            c = comp[rd]
            r = dict(kind=kind, readout=rd, cohorts=len(units), patients=n, null_mean=float(nul.mean()),
                     null_sd=float(nul.std()), null_q95=float(np.quantile(nul, .95)), mde=float(MDE_K * nul.std()),
                     decidable_mde=bool(MDE_K * nul.std() <= MAX_MDE))
            for f, i in (("matched", 0), ("responsive", 1), ("responsive_broad", 2)):
                s = c[i * N_SET:(i + 1) * N_SET]
                r[f"{f}_median"], r[f"{f}_q95"] = float(np.nanmedian(s)), float(np.nanquantile(s, .95))
            for rho in (0.3, 0.5):
                r[f"phase_halfwidth_h_at_{rho}"] = float(phase_halfwidth_h(rho, n))
            rows.append(r)
            log(f"part 1 {kind} {rd}: {n} patients, null {r['null_mean']:.3f} +- {r['null_sd']:.3f}")
    return pd.DataFrame(rows)


# ---- Part 2: full model -----------------------------------------------------------------------------------------

def cohort_table(study):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        pb = M.table(study)
    except Exception as e:  # noqa: BLE001  (a cohort that cannot be built is reported, as A346)
        print(f"{study} not built: {type(e).__name__}: {e}", flush=True)
        return None
    keep = set(ANCHORS) | set(TARGETS.values()) | set(NORMAL.values())
    pb = pb[pb.compartment.isin(keep)].copy()
    pb["donor"] = study + "|" + pb.donor.astype(str)
    pb["stratum"] = study + "|" + (pb.stratum.astype(str) if "stratum" in pb else "0")
    return pb


def pooled_table(jobs, smoke):
    from joblib import Parallel, delayed
    cohorts = M.cohort_list(smoke)
    T = [t for t in Parallel(n_jobs=jobs)(delayed(cohort_table)(s) for s in cohorts) if t is not None and len(t)]
    return pd.concat(T, ignore_index=True)


def part2_power(table, smoke):
    rows, gates, design = [], [], []
    for rd, wf in READOUTS.items():
        ds = cs.Dataset(table, anchors=ANCHORS, waveform=wf, mixed_strata=True, name=f"A387_{rd}")
        g = ds.subset(compartments=ANCHORS).gate(n_perm=100 if smoke else 1000, seed=387)
        gates.append(dict(readout=rd, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share,
                          share_null_mean=g.share_null_mean, n_pairs=g.n_pairs, n_donors=g.n_donors))
        log(f"gate {rd}: p_pairing {g.p_pairing:.4f}, p_acrophase {g.p_acrophase:.4f}")
        hosted = set(ds.host().index.astype(str))
        for kind, k in TARGETS.items():
            donors_k = set(ds.table.resection[ds.table.compartment == k].astype(str))
            design.append(dict(readout=rd, kind=kind, target=k, patients=len(donors_k),
                               patients_with_host=len(donors_k & hosted),
                               cohorts=len({d.split("|")[0] for d in donors_k})))
            for za in (1.0, 0.5):
                sim = cs.simulate.like(ds, k, loading=za)
                rng = np.random.default_rng([387, KCODE[kind], int(za * 10), rd == "core"])
                reps = {c: [sim.draw(c, rng) for _ in range(N_SIM)] for c in C_GRID}
                thr = float(np.nanquantile([r["coupling"] for r in reps[C_GRID[0]]], .95))
                for c in C_GRID:
                    cc = np.array([r["coupling"] for r in reps[c]])
                    ph = np.array([np.angle(r["z"]) for r in reps[c]])
                    sd_h = float(np.sqrt(-2 * np.log(np.abs(np.nanmean(np.exp(1j * ph)))))) * HOURS
                    rows.append(dict(readout=rd, kind=kind, z_abs=za, c_true=c, c_median=float(np.nanmedian(cc)),
                                     null_threshold=thr, power=float(np.nanmean(cc > thr)),
                                     offset_halfwidth90_h=1.645 * sd_h))
            log(f"part 2 {kind} {rd}: simulated")
    P = pd.DataFrame(rows)
    mde = []
    for (rd, kind, za), g in P.groupby(["readout", "kind", "z_abs"]):
        ok = g[(g.c_true > C_GRID[0]) & (g.power >= 0.8)]
        h = g[g.c_true == 0.5].offset_halfwidth90_h.iloc[0]
        mde.append(dict(readout=rd, kind=kind, z_abs=za, mde_c=float(ok.c_true.min()) if len(ok) else np.nan,
                        offset_halfwidth90_h_at_c05=h,
                        decidable_coupling=bool(len(ok) and ok.c_true.min() <= MAX_MDE),
                        decidable_offset=bool(h <= BOUND_H)))
    return P, pd.DataFrame(mde), pd.DataFrame(gates), pd.DataFrame(design)


def stage_power(args):
    out = OUT / ("power_smoke" if args.smoke else "power")
    out.mkdir(parents=True, exist_ok=True)
    p1 = part1_power(args.smoke)
    p1.to_csv(out / "part1_power.csv", index=False)
    log("part 1 done")
    table = pooled_table(args.jobs, args.smoke)
    counts = table.groupby("compartment").donor.nunique().rename("patients").reset_index()
    counts.to_csv(out / "part2_counts.csv", index=False)
    P, mde, gates, design = part2_power(table, args.smoke)
    P.to_csv(out / "part2_power_curve.csv", index=False)
    mde.to_csv(out / "part2_power.csv", index=False)
    gates.to_csv(out / "part2_gate.csv", index=False)
    design.to_csv(out / "part2_design.csv", index=False)
    pd.set_option("display.width", 250)
    for name, t in (("part 1", p1), ("counts", counts), ("gate", gates), ("design", design), ("part 2", mde)):
        log(f"{name}\n" + t.round(4).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power",), required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    global N_PERM, N_DRAW, N_SIM
    if args.smoke:
        N_PERM, N_DRAW, N_SIM = 200, 20, 20
    t0 = time.time()
    stage_power(args)
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
