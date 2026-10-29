"""A155: bias from reading anchors on a host that contains them (registered in paper/docs/registered_experiments.md, A155).

`Dataset.fit` reads every compartment's relative amplitude |z| and offset arg z on the full anchor host; an anchor's own noise
is then part of its regressor. Here, in simulated cohorts with known amplitude, offset and coupling, each anchor is read on
(a) the full host (as `Dataset.fit`: `host_a` of all anchors, `on_host`) and (b) a leave-one-out host (the host fitted on the
other anchors, then put on the full host's gauge; see README). Non-anchors are read on the full host in both (reference).

Designs (cells): m anchors in {3, 4, 6, 10} + 4 non-anchor targets, n donors in {30, 100}, depth single-cell or bulk
(16 cells), and the A123-like design (40 compartments, all anchors, n = 111, bulk), 200 replicates, seed 155.
Truth from A123 truth (i) (amplitude, offset) and the A123 per-tissue coupling (`organ_map.csv`); see README.

    RUN_CAP_GB=10 paper/common/run.sh paper/exploratory/A155_anchor_self_inclusion/scripts/compute_a155.py --jobs 6
    options: --reps N (default 200), --out DIR (default ../output), --cells name,name (default all)
Outputs (DIR): raw/<cell>.csv.gz (one row per replicate x compartment: truth, estimates under a and b), summary_bias.csv,
summary_a123like.csv, decision.csv. A rerun resumes from the raw files (a cell is recomputed only if its file is missing).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

import depth_values  # noqa: E402

OUT = HERE.parent / "output"
A123 = PAPER / "exploratory" / "A123_working_hours" / "output"
HOURS = 24 / (2 * np.pi)
SEED = 155
WINDOW, JITTER = (8.0, 18.0), 1.0  # A17 / A123 sampling spread (the simulator's default window)
M_GRID, N_GRID, DEPTHS, N_TARGETS = (3, 4, 6, 10), (30, 100), ("single_cell", "bulk"), 4
A123_LIKE = "a123like"


def cells() -> dict:
    out = {f"m{m}_n{n}_{dp}": dict(m=m, n=n, depth=dp, targets=N_TARGETS)
           for m in M_GRID for n in N_GRID for dp in DEPTHS}
    out[A123_LIKE] = dict(m=40, n=111, depth="bulk", targets=0)
    return out


# ---- truth ----------------------------------------------------------------------------------------------------------

def truth_table() -> pd.DataFrame:
    """The 40 A123 tissues: amplitude and offset (h) of truth (i), coupling of the A123 estimate (clipped to [0.15, 0.95])."""
    t = pd.read_csv(A123 / "truth_i.csv").set_index("compartment")[["amplitude", "offset_h"]]
    o = pd.read_csv(A123 / "organ_map.csv").set_index("tissue").coupling.clip(0.15, 0.95)
    return t.join(o.rename("coupling")).sort_index()


def psi_for(coupling: float, z_abs: float) -> float:
    """Clock-plane noise psi giving true coupling `coupling` for a loading of modulus `z_abs` at the sampling spread
    (psi = |z|^2 V (1 - c) / (2 c), V = trace of the covariance of (cos, sin) of the sampled time)."""
    t = np.random.default_rng(0).uniform(*WINDOW, 200000) + np.random.default_rng(1).normal(0, JITTER, 200000)
    V = np.trace(np.cov(np.column_stack([np.cos(t / HOURS), np.sin(t / HOURS)]).T))
    return z_abs ** 2 * V * (1 - coupling) / (2 * coupling)


def design(cell: str) -> pd.DataFrame:
    """Compartments of a cell: role, true loading z, true coupling, psi. Anchors a00.., targets t0..; the 40 tissues are
    permuted once (seed 155): the first `targets` are targets, the next m anchors (nested across m)."""
    p = cells()[cell]
    T = truth_table()
    if cell == A123_LIKE:
        order, names = T.index.to_list(), [f"a{i:02d}" for i in range(40)]
        roles = ["anchor"] * 40
    else:
        order = list(np.random.default_rng(SEED).permutation(T.index.to_numpy()))
        tg, an = order[:p["targets"]], order[p["targets"]:p["targets"] + p["m"]]
        order = an + tg
        names = [f"a{i:02d}" for i in range(p["m"])] + [f"t{i}" for i in range(p["targets"])]
        roles = ["anchor"] * p["m"] + ["target"] * p["targets"]
    d = T.loc[order].reset_index(drop=True)
    d.insert(0, "compartment", names)
    d["role"] = roles
    d["tissue"] = order
    d["z"] = d.amplitude * np.exp(1j * d.offset_h / HOURS)
    d["psi"] = [psi_for(c, a) for c, a in zip(d.coupling, d.amplitude)]
    return d


# ---- one replicate --------------------------------------------------------------------------------------------------

def centred(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Log amplitude minus its mean and offset (h) minus the circular mean over compartments: the reference-free
    readout of `Dataset.fit`."""
    la = np.log(np.abs(z))
    ph = np.angle(z)
    mean = np.angle(np.exp(1j * ph).mean())
    off = (np.angle(np.exp(1j * (ph - mean)))) * HOURS
    return la - la.mean(), off


def read_hosts(d: A.Data, anchors: list, comps: list) -> dict:
    """z and coupling of each compartment on (a) the full anchor host and (b) the leave-one-out host (anchors) put on the
    gauge of the full host; non-anchors: full host for both."""
    full = A.host(d, anchors)
    ref = full.ref
    out = {}
    for k in comps:
        ra = A.on_host(d, k, full)
        rb = ra
        if k in anchors:
            h = A.host(d, anchors, exclude=k, ref=ref if k != ref else None)
            rb = A.on_host(d, k, h)
            if h.ref != ref:  # the reference compartment left out: rescale to the gauge of the full host
                zl = np.array([h.z[list(h.compartments).index(j)] for j in h.compartments])
                zf = np.array([full.z[list(full.compartments).index(j)] for j in h.compartments])
                s = np.sum(np.conj(zl) * zf) / np.sum(np.abs(zl) ** 2)
                rb = dict(rb, z=rb["z"] * s)
        out[k] = (ra, rb)
    return out


def replicate(cell: str, rep: int) -> pd.DataFrame:
    p, D = cells()[cell], design(cell)
    dv = depth_values.load()
    seed = int(np.random.SeedSequence([SEED, list(cells()).index(cell), rep]).generate_state(1)[0])
    comps = {r.compartment: Compartment(z=complex(r.z), psi=r.psi) for r in D.itertuples()}
    tab, _ = cs.simulate(cs.gtex(), comps, n=p["n"], window_h=WINDOW, jitter_h=JITTER,
                         depth=depth_values.depth_for(dv[p["depth"]]), seed=seed)
    anchors = list(D.compartment[D.role == "anchor"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = cs.Dataset(tab.rename(columns={"resection": "donor"}), anchors=anchors)
    est = read_hosts(ds.data, ds.anchors, list(D.compartment))
    zt = D.z.to_numpy()
    lat, ot = centred(zt)
    rows = dict(cell=cell, rep=rep, compartment=D.compartment, role=D.role, true_log_amp=lat, true_offset_h=ot,
                true_coupling=D.coupling)
    for j, tag in enumerate(("a", "b")):
        z = np.array([est[k][j]["z"] for k in D.compartment])
        la, off = centred(z)
        rows[f"log_amp_{tag}"], rows[f"offset_h_{tag}"] = la, off
        rows[f"coupling_{tag}"] = [est[k][j]["coupling"] for k in D.compartment]
    return pd.DataFrame(rows)


def chunk(cell: str, reps: list) -> pd.DataFrame:
    limit_threads(1)
    return pd.concat([replicate(cell, r) for r in reps])


# ---- scoring --------------------------------------------------------------------------------------------------------

def wrap12(x):
    return (np.asarray(x) + 12.0) % 24.0 - 12.0


def circ_corr(a, b) -> float:
    """Jammalamadaka-SenGupta circular correlation of two angle vectors (h)."""
    a, b = np.asarray(a) / HOURS, np.asarray(b) / HOURS
    sa, sb = np.sin(a - np.angle(np.exp(1j * a).sum())), np.sin(b - np.angle(np.exp(1j * b).sum()))
    return float((sa * sb).sum() / np.sqrt((sa ** 2).sum() * (sb ** 2).sum()))


def score(out: Path):
    raw = pd.concat([pd.read_csv(f) for f in sorted((out / "raw").glob("*.csv.gz"))])
    rows = []
    for (cell, role), g in raw.groupby(["cell", "role"]):
        for tag in ("a", "b"):
            e_la = g[f"log_amp_{tag}"] - g.true_log_amp
            e_off = wrap12(g[f"offset_h_{tag}"] - g.true_offset_h)
            e_c = g[f"coupling_{tag}"] - g.true_coupling
            per = pd.DataFrame({"comp": g.compartment, "la": e_la, "off": e_off, "c": e_c}).groupby("comp")
            bias = per.mean()
            reps = g.rep.nunique()
            rows.append(dict(cell=cell, role=role, host=tag, reps=reps, n_comp=len(bias),
                             amp_bias_pct_signed=100 * np.expm1(bias.la.mean()),
                             amp_bias_pct_mean_abs=100 * np.mean(np.abs(np.expm1(bias.la))),
                             amp_bias_pct_max_abs=100 * np.max(np.abs(np.expm1(bias.la))),
                             amp_bias_pct_se=100 * float(per.la.std().mean() / np.sqrt(reps)),
                             offset_bias_h_signed=bias.off.mean(), offset_bias_h_mean_abs=np.mean(np.abs(bias.off)),
                             offset_bias_h_max_abs=np.max(np.abs(bias.off)),
                             offset_bias_h_se=float(per.off.std().mean() / np.sqrt(reps)),
                             coupling_bias_signed=bias.c.mean(), coupling_bias_mean_abs=np.mean(np.abs(bias.c)),
                             rmse_log_amp=float(np.sqrt((e_la ** 2).mean())), rmse_offset_h=float(np.sqrt((e_off ** 2).mean())),
                             rmse_coupling=float(np.sqrt((e_c ** 2).mean()))))
    S = pd.DataFrame(rows)
    S.to_csv(out / "summary_bias.csv", index=False)
    # A123-like: circular correlation of offsets and Pearson of log2 amplitudes against truth, per replicate
    g = raw[raw.cell == A123_LIKE]
    if len(g):
        r = []
        for rep, x in g.groupby("rep"):
            row = dict(rep=rep)
            for tag in ("a", "b"):
                row[f"offset_r_{tag}"] = circ_corr(x[f"offset_h_{tag}"], x.true_offset_h)
                row[f"amp_r_{tag}"] = float(np.corrcoef(x[f"log_amp_{tag}"], x.true_log_amp)[0, 1])
            r.append(row)
        R = pd.DataFrame(r)
        R.to_csv(out / "a123like_per_rep.csv", index=False)
        d_off, d_amp = R.offset_r_b - R.offset_r_a, R.amp_r_b - R.amp_r_a
        pd.DataFrame([dict(reps=len(R), offset_r_a=R.offset_r_a.mean(), offset_r_b=R.offset_r_b.mean(),
                           offset_r_change_b_minus_a=d_off.mean(), offset_r_change_se=d_off.std() / np.sqrt(len(R)),
                           amp_r_a=R.amp_r_a.mean(), amp_r_b=R.amp_r_b.mean(),
                           amp_r_change_b_minus_a=d_amp.mean(), amp_r_change_se=d_amp.std() / np.sqrt(len(R)))]
                     ).to_csv(out / "summary_a123like.csv", index=False)
    # decision rule: anchor bias > 5 % in |z| or > 0.25 h in offset under (a), A17-typical designs (m 3-6, n 30-100)
    a = S[(S.host == "a") & (S.role == "anchor") & (S.cell != A123_LIKE)].copy()
    a["m"] = a.cell.str.extract(r"m(\d+)_").astype(int)
    a["n"] = a.cell.str.extract(r"_n(\d+)_").astype(int)
    a["depth"] = a.cell.str.extract(r"_n\d+_(.*)")
    a["typical"] = a.m.between(3, 6)
    a["exceeds_amp_5pct"] = a.amp_bias_pct_mean_abs > 5
    a["exceeds_offset_0.25h"] = a.offset_bias_h_mean_abs > 0.25
    a["exceeds_amp_5pct_signed"] = a.amp_bias_pct_signed.abs() > 5
    a["exceeds_offset_0.25h_signed"] = a.offset_bias_h_signed.abs() > 0.25
    a[["cell", "m", "n", "depth", "typical", "amp_bias_pct_signed", "amp_bias_pct_mean_abs", "offset_bias_h_signed",
       "offset_bias_h_mean_abs", "exceeds_amp_5pct", "exceeds_offset_0.25h", "exceeds_amp_5pct_signed",
       "exceeds_offset_0.25h_signed"]].to_csv(out / "decision.csv", index=False)


def main():
    from common import cli
    ap = cli.parser(__doc__)
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--cells", default=None)
    ap.add_argument("--block", type=int, default=10, help="replicates per task")
    args = ap.parse_args()
    out = args.out
    (out / "raw").mkdir(parents=True, exist_ok=True)
    todo = [c for c in (args.cells.split(",") if args.cells else cells()) if not (out / "raw" / f"{c}.csv.gz").exists()]
    # big cells first (the 40-compartment design dominates)
    todo.sort(key=lambda c: -cells()[c]["m"])
    tasks = [(c, list(range(i, min(i + args.block, args.reps)))) for c in todo for i in range(0, args.reps, args.block)]
    print(f"{len(todo)} cells, {len(tasks)} tasks, {args.jobs} jobs", flush=True)
    parts, t0, last = {c: [] for c in todo}, time.time(), time.time()
    if tasks:
        with cf.ProcessPoolExecutor(args.jobs) as ex:
            futs = {ex.submit(chunk, c, r): (c, r) for c, r in tasks}
            for f in cf.as_completed(futs):
                c, r = futs[f]
                parts[c].append(f.result())
                if len(parts[c]) == -(-args.reps // args.block):
                    pd.concat(parts[c]).sort_values(["rep", "compartment"]).to_csv(out / "raw" / f"{c}.csv.gz", index=False)
                    parts[c] = []
                    print(f"{time.strftime('%H:%M:%S')} cell {c} done ({time.time() - t0:.0f} s)", flush=True)
                elif time.time() - last > 600:
                    print(f"{time.strftime('%H:%M:%S')} {sum(len(v) for v in parts.values())} blocks pending cells", flush=True)
                    last = time.time()
    score(out)
    print(f"done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
