"""A139 library: simulated design, E-a (clocksync `Dataset.gene_loadings`, E4) and E-b (iterative leave-one-gene-out
waveform refit) with donor-bootstrap flags. Registered in paper/docs/registered_experiments.md (A139)."""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT))

import clocksync as cs  # noqa: E402
from clocksync._core.plane import log_fraction  # noqa: E402
from clocksync._core.study import Study  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parent / "output"
HOURS = 24.0 / (2 * np.pi)
WF = cs.gtex()
GENES = list(WF.genes)
ANCHORS = ["A1", "A2", "A3", "A4"]
TARGET = "T"
C_ANCHOR = {"A1": 0.8, "A2": 0.9, "A3": 0.9, "A4": 1.0}
C_TARGET, Z_TARGET = 0.8, 1.0 * np.exp(1j * 3.0 / HOURS)  # target: coupling 0.8, amplitude 1, offset +3 h (fixed choice)
SHIFT_H = 6.0
SHIFTED = {"S0": [], "S1": ["PER2"], "S3": ["PER2", "TEF", "NR1D2"], "SW": ["PER2", "TEF", "NR1D2"]}
WINDOWS = {"24h": (0.0, 24.0), "08-18h": (8.0, 18.0), "10-16h": (10.0, 16.0)}
JITTER = 1.0
NS = (40, 80, 200)
SCENARIOS = ("S0", "S1", "S3", "SW")
N_BOOT, N_PERM_LR = 200, 500


def host_variance(window, jitter=JITTER) -> float:
    r = np.random.default_rng(0)
    t = r.uniform(*window, 400_000) + (r.normal(0, jitter, 400_000) if window != WINDOWS["24h"] else 0.0)
    return float(np.trace(np.cov(np.column_stack([np.cos(t / HOURS), np.sin(t / HOURS)]).T)))


def psi_of(c, z, V):
    return z * z * V * (1 - c) / (2 * c)


def cohort(window_key, n, scen, rep):
    """Long table (`donor`, `compartment`, `depth`, genes). Same seed for the four scenarios of a (window, n, rep)."""
    win = WINDOWS[window_key]
    V = host_variance(win)
    sh = {g: SHIFT_H for g in SHIFTED[scen]}
    comps = {k: Compartment(z=1.0, psi=psi_of(c, 1.0, V), gene_shift_h=sh if scen == "SW" else None)
             for k, c in C_ANCHOR.items()}
    comps[TARGET] = Compartment(z=Z_TARGET, psi=psi_of(C_TARGET, abs(Z_TARGET), V),
                                gene_shift_h=sh if scen != "S0" else None)
    seed = int(np.random.SeedSequence([139, list(WINDOWS).index(window_key), n, rep]).generate_state(1)[0])
    tab, _ = cs.simulate(WF, comps, n=n, window_h=win, jitter_h=JITTER, seed=seed)
    return tab.rename(columns={"resection": "donor"})


def truth_dev(scen, centred=False):
    """Nominal true deviation (h) of every gene from the conserved waveform in the target; `centred`: deviation of the
    gene-centred loading (the quantity E4 reads)."""
    d = np.array([SHIFT_H if g in SHIFTED[scen] else 0.0 for g in GENES])
    if not centred:
        return d
    X = WF.amp * np.exp(1j * (WF.psi + d / HOURS))
    X0 = WF.amp * np.exp(1j * WF.psi)
    return np.angle((X - X.mean()) / (X0 - X0.mean())) * HOURS


def wrap(h):
    return (np.asarray(h, float) + 12.0) % 24.0 - 12.0


def resample_table(tab, rng):
    donors = np.array(sorted(tab.donor.unique()))
    pick = rng.choice(donors, len(donors))
    parts = []
    for j, d in enumerate(pick):
        parts.append(tab[tab.donor == d].assign(donor=f"{d}#{j}"))
    return pd.concat(parts, ignore_index=True)


# ---- E-a --------------------------------------------------------------------------------------------------------

def e_a(tab, seed, with_gate=True):
    """E4 on the target with the anchors as host: dev_h per gene, LR p, and donor-bootstrap interval of dev_h."""
    ds = cs.Dataset(tab, anchors=ANCHORS, name="a139")
    genes, summ = ds.gene_loadings(TARGET, n_perm=N_PERM_LR, seed=seed, hosts=ANCHORS)
    rng = np.random.default_rng([139, seed])
    B = []
    for _ in range(N_BOOT):
        try:
            g, _s = cs.Dataset(resample_table(tab, rng), anchors=ANCHORS).gene_loadings(TARGET, n_perm=0, hosts=ANCHORS)
            B.append(g.dev_h.to_numpy())
        except Exception:  # a degenerate resample is skipped (counted in n_boot_ok)
            continue
    B = np.array(B)
    d0 = genes.dev_h.to_numpy()
    Bw = d0 + wrap(B - d0) if len(B) else np.full((1, len(d0)), np.nan)
    lo, hi = np.nanquantile(Bw, 0.025, 0), np.nanquantile(Bw, 0.975, 0)
    fit = ds.fit(n_boot=0, seed=seed, force=True)
    gate = ds.gate(n_perm=200, seed=seed) if with_gate else None
    t = fit.table.loc[TARGET]
    return dict(dev=d0, se=genes.dev_se_h.to_numpy(), lo=lo, hi=hi, lr_p=float(summ.p), n_boot_ok=len(B),
                z_amp=float(summ.z_amplitude), z_off_h=float(summ.z_offset_h),
                fit_logamp=float(t.log_amp_rel), fit_off=float(t.offset_h),
                gate_share=float(gate.clock_share) if gate else np.nan, gate_pass=bool(gate.passes) if gate else None,
                gate_p=float(gate.p_acrophase) if gate else np.nan)


# ---- E-b --------------------------------------------------------------------------------------------------------

def _wf(X, genes, name="refit"):
    return Waveform(tuple(genes), np.abs(X), np.angle(X), None, None, name)


def _pool_loading(tab, X, g, comps, all_comps):
    """Loading of gene g (uncentred, complex) from the host of the other 13 genes over all compartments, pooled over
    compartments: per compartment y_kg regressed on the host state (one intercept), c_k = conj(z_k) Xc_g, combined by
    inverse-variance complex least squares; then X_g = (14 Xc_g + sum of the other 13 X) / 13."""
    others = [j for j in range(len(GENES)) if j != g]
    wf13 = _wf(X[others], [GENES[j] for j in others])
    st = Study(tab.rename(columns={"donor": "resection"}), wf13, ref="A1")
    ft = st.fit()
    pos = {r: i for i, r in enumerate(ft.resections)}
    y = log_fraction(tab[GENES].to_numpy(float), tab.depth.to_numpy(float))
    y -= y.mean(1, keepdims=True)
    num = den = 0.0
    for kk, k in enumerate(ft.compartments):
        i = np.flatnonzero(tab.compartment.to_numpy() == k)
        f = ft.f[[pos[r] for r in tab.donor.to_numpy()[i]]]
        D = np.column_stack([np.ones(len(i)), f])
        th, *_ = np.linalg.lstsq(D, y[i, g], rcond=None)
        s2 = np.sum((y[i, g] - D @ th) ** 2) / max(len(i) - 3, 1)
        v = s2 * np.linalg.inv(D.T @ D)[1:, 1:].trace() / 2.0
        c = th[1] + 1j * th[2]
        zk = ft.z[kk]
        num += zk * c / v
        den += abs(zk) ** 2 / v
    Xc = num / den
    return (len(GENES) * Xc + X[others].sum()) / (len(GENES) - 1)


def refit(tab, X0, tol_h=0.1, max_rounds=10):
    """Iterative leave-one-gene-out refit: in each round every gene is refitted against the host of the other 13 genes
    at their current values (simultaneous update; the sequential variant drifts in scale and was not used), then the new
    waveform is divided by its complex least-squares scale on the previous one (the overall amplitude and rotation are
    not identified by the refit and otherwise drift). Stops when
    no acrophase changes by >= tol_h. Returns X, rounds, converged."""
    X = X0.copy()
    for r in range(1, max_rounds + 1):
        new = np.array([_pool_loading(tab, X, g, None, None) for g in range(len(GENES))])
        new = new / (np.vdot(X, new) / np.vdot(X, X))  # gauge: complex scale (amplitude, rotation) of the whole waveform
        change = np.abs(wrap(np.angle(new / X) * HOURS)).max()
        X = new
        if change < tol_h:
            return X, r, True
    return X, max_rounds, False


def e_b(tab, seed, boot_rounds=3):
    """Waveform refit; deviation = arg(X_refit / X_gtex) (h). Bootstrap replicates are warm-started at the full-data
    refit (same stopping rule)."""
    X0 = WF.amp * np.exp(1j * WF.psi)
    X, rounds, conv = refit(tab, X0)
    d0 = np.angle(X / X0) * HOURS
    rng = np.random.default_rng([139, 1, seed])
    B = []
    for _ in range(N_BOOT):
        try:
            Xb, _, _ = refit(resample_table(tab, rng), X, max_rounds=boot_rounds)
            B.append(np.angle(Xb / X0) * HOURS)
        except Exception:
            continue
    B = np.array(B)
    Bw = d0 + wrap(B - d0) if len(B) else np.full((1, len(d0)), np.nan)
    return dict(dev=wrap(d0), lo=np.nanquantile(Bw, 0.025, 0), hi=np.nanquantile(Bw, 0.975, 0), rounds=rounds,
                converged=conv, n_boot_ok=len(B))


def run_rep(window_key, n, scen, rep):
    """Gene-level rows and one replicate-level row of a (window, n, scenario, replicate)."""
    tab = cohort(window_key, n, scen, rep)
    seed = int(np.random.SeedSequence([139, 2, list(WINDOWS).index(window_key), n, rep]).generate_state(1)[0])
    base = dict(window=window_key, n=n, scenario=scen, rep=rep)
    a = e_a(tab, seed, with_gate=scen in ("S3", "SW"))
    rows = pd.DataFrame(dict(**base, gene=GENES, shifted=[g in SHIFTED[scen] for g in GENES],
                             truth_h=truth_dev(scen), truth_centred_h=truth_dev(scen, True),
                             a_dev=a["dev"], a_se=a["se"], a_lo=a["lo"], a_hi=a["hi"]))
    rows["a_flag_boot"] = (rows.a_lo > 0) | (rows.a_hi < 0)
    rows["a_flag_lr"] = a["lr_p"] < 0.05
    summ = dict(**base, a_lr_p=a["lr_p"], a_n_boot_ok=a["n_boot_ok"], a_z_amp=a["z_amp"], a_z_off_h=a["z_off_h"],
                fit_logamp=a["fit_logamp"], fit_off_h=a["fit_off"], gate_share=a["gate_share"], gate_pass=a["gate_pass"],
                gate_p=a["gate_p"])
    if scen in ("S0", "SW"):
        b = e_b(tab, seed)
        rows["b_dev"], rows["b_lo"], rows["b_hi"] = b["dev"], b["lo"], b["hi"]
        rows["b_flag_boot"] = (rows.b_lo > 0) | (rows.b_hi < 0)
        summ.update(b_rounds=b["rounds"], b_converged=b["converged"], b_n_boot_ok=b["n_boot_ok"])
    return rows, summ
