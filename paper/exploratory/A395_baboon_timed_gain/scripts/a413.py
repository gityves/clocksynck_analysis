"""A413 (do an organ's rhythmic outputs inherit its own clock's departure from clock time?), registered 2026-11-27.

Baboon (Mure 2018, 12 animals, ZT known). Per sample s (animal i, organ k) of a gene set read with its own baboon-fitted
waveform: r_s = z_s - zhat_s from `Dataset.timed_loadings("zt", return_residuals=True)`, relative departure
delta_s = r_s / zhat_s. Local part = delta_s minus the animal's mean delta over its other passing organs (leave-one-out),
divided by the organ's RMS. rho_j = Re complex correlation of the clock's and output set j's local parts over jointly
passing samples; rho_med = median over the A407 output sets with >= 15 jointly passing organs. Companions: 200
expression-matched non-rhythmic 14-gene pseudo-clocks (A263 rule, seed [413, 0]) read at the clock positions with the
clock's waveform; their residual is normalised by the clock's zhat_s (clock's passing organs and samples).

Stages:
  power  precondition counts, companion null (companions.csv: contains no clock-output value), injection power
         (output residual = zhat_j x dl_c, clock local departure, + its own residuals permuted between animals within
         organ; 100 draws; full and half inheritance; zero inheritance as a calibration check). The real clock-output
         correlation is not computed here.
  score  H413.1 (rho_med of the clock vs companion q95 and 1,000 permutations of the clock residuals between animals
         within organ) and the descriptive readouts. NOT RUN before the power gate is read.

    RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A395_baboon_timed_gain/scripts/a413.py --stage power --jobs 4
Outputs (`../output/a413/`): outputs.npz (per output set: residual and zhat grids of its passing organs), companions.csv,
power.csv, power_injection.csv; tests.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a407 as M  # noqa: E402  (pool, sets, set_waveform, perms, load_genome)
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.timed import W  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

B, A = M.B, M.A
OUT = HERE.parents[1] / "output" / "a413"
SEED = 413
KEYS, WF, TRIO = M.KEYS, M.WF, M.TRIO
N_SETS, N_COMP, N_DRAW, N_PERM, MIN_JOINT = 200, 200, 100, 1000, 15
S = {}


# ---- grids: animals x organs, complex, NaN where absent --------------------------------------------------------

def grids(ds, zt, per, animals, organs=None):
    """Residual r and fitted zhat of every sample as (animal, organ) grids; organs None: all those fitted."""
    tab, res = B.timed(ds, zt, residuals=True)
    organs = list(tab.index) if organs is None else organs
    a_k = tab.abs_amplitude * np.exp(1j * W * tab.clock_offset_h)
    r = res.reset_index()
    r.columns = ["donor", "compartment", "res"]
    r["donor"] = r.donor.astype(str).str.split("#").str[0]
    r = r[r.compartment.isin(organs)]
    ai = pd.Index(animals).get_indexer(r.donor)
    ki = pd.Index(per).get_indexer(r.compartment)
    R = np.full((len(animals), len(per)), np.nan + 0j)
    Z = np.full_like(R, np.nan)
    R[ai, ki] = r.res.to_numpy()
    Z[ai, ki] = (a_k.reindex(r.compartment).to_numpy() * np.exp(1j * W * zt.reindex(r.donor).to_numpy(float)))
    return R, Z


def passing(ds, zt, per):
    obs = B.timed(ds, zt).abs_amplitude.reindex(per)
    null = np.array([B.timed(ds, p).abs_amplitude.reindex(per).to_numpy(float) for p in S["perms"]])
    return obs.to_numpy(float) > np.nanquantile(null, 0.95, axis=0)


def local(D):
    """Leave-one-out local part (animal's mean over its other organs removed) of a grid; single-organ animals NaN."""
    ok = np.isfinite(D)
    n = ok.sum(1, keepdims=True)
    tot = np.where(ok, D, 0).sum(1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = D - (tot - np.where(ok, D, 0)) / (n - 1)
    out[np.broadcast_to(n < 2, D.shape)] = np.nan
    return out


def standardise(L):
    """Divide each organ column by its RMS."""
    with np.errstate(invalid="ignore"):
        rms = np.sqrt(np.nanmean(np.abs(L) ** 2, axis=0, keepdims=True))
    return L / rms


def rho(a, b, part="both"):
    """Re complex correlation sum a conj(b) / sqrt(sum|a|^2 sum|b|^2) over jointly finite entries (part: Re / Im only)."""
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return np.nan
    x, y = a[m], b[m]
    if part == "re":
        x, y = x.real + 0j, y.real + 0j
    elif part == "im":
        x, y = 1j * x.imag, 1j * y.imag
    return float((x * np.conj(y)).sum().real / np.sqrt((np.abs(x) ** 2).sum() * (np.abs(y) ** 2).sum()))


def rho_med(lc, outs):
    """rho of local part lc against each output set's standardised local part; median over sets with >= 15 joint organs."""
    v = []
    for lj in outs:
        if (np.isfinite(lc) & np.isfinite(lj)).any(0).sum() >= MIN_JOINT:
            v.append(rho(lc, lj))
    return float(np.nanmedian(v)) if v else np.nan, len(v)


# ---- setup ------------------------------------------------------------------------------------------------------

def setup():
    pb, per, zt, genes = M.load_genome()
    animals = sorted(zt.index)
    clock = list(WF.genes)
    wf_c = M.set_waveform(pb[KEYS + clock], per, zt, clock)
    S.update(pb=pb, per=per, zt=zt, animals=animals, wf_c=wf_c, perms=M.perms(zt))
    ds = cs.Dataset(pb[KEYS + clock], anchors=per, waveform=wf_c, name="a413 clock", min_donors=5)
    S["ds_c"] = ds
    S["pass_c"] = passing(ds, zt, per)
    organs = [k for k, p in zip(per, S["pass_c"]) if p]
    S["Rc"], S["Zc"] = grids(ds, zt, per, animals, organs)
    return organs


def output_task(i):
    t0 = time.time()
    s = S["sets"][S["sets"].set == i].sort_values("position")
    genes = list(s.gene)
    pb, per, zt = S["pb"][KEYS + genes], S["per"], S["zt"]
    wf = M.set_waveform(pb, per, zt, genes)
    ds = cs.Dataset(pb, anchors=per, waveform=wf, name=f"a413 set {i}", min_donors=5)
    p = passing(ds, zt, per)
    R, Z = grids(ds, zt, per, S["animals"], [k for k, q in zip(per, p) if q])
    return i, p, R, Z, time.time() - t0


def companion_task(i):
    genes = list(S["comp"][i].genes)
    wf = Waveform(tuple(genes), np.asarray(S["wf_c"].amp), np.asarray(S["wf_c"].psi), None, None, f"pseudo{i}")  # as_set
    pb = S["pbc"][KEYS + genes]
    ds = cs.Dataset(pb, anchors=S["per"], waveform=wf, name=f"a413 pseudo {i}", min_donors=5)
    organs = [k for k, p in zip(S["per"], S["pass_c"]) if p]
    R, _ = grids(ds, S["zt"], S["per"], S["animals"], organs)
    return i, genes, R


def output_sets(jobs):
    f = OUT / "outputs.npz"
    if f.exists():
        z = np.load(f)
        return z["P"], z["R"], z["Z"]
    t0 = time.time()
    sets = pd.read_csv(M.OUT / "sets.csv")
    S["sets"] = sets
    P = np.zeros((N_SETS, len(S["per"])), bool)
    R = np.full((N_SETS, len(S["animals"]), len(S["per"])), np.nan + 0j)
    Z = np.full_like(R, np.nan)
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, (i, p, r, z, sec) in enumerate(pool.imap_unordered(output_task, range(N_SETS))):
            P[i], R[i], Z[i] = p, r, z
            if (k + 1) % 5 == 0 or k + 1 == N_SETS:
                print(f"outputs {k + 1}/{N_SETS} {time.time() - t0:.0f}s (last {sec:.0f}s)", flush=True)
    np.savez_compressed(f, P=P, R=R, Z=Z, organs=np.array(S["per"]), animals=np.array(S["animals"]))
    return P, R, Z


def out_locals(P, R, Z):
    """Standardised local parts of the output sets."""
    return [standardise(local(R[j] / Z[j])) for j in range(len(R))]


def companions(jobs):
    f = OUT / "companions.csv"
    lm = A.log_mean("baboon")
    sets = pseudo_waveforms(lm, N_COMP, np.random.default_rng([SEED, 0]), WF, "matched")  # A263 rule, baboon table
    S["comp"] = sets
    genes = sorted({g for w in sets for g in w.genes})
    pbc = A.table("baboon", genes)
    pbc["donor"] = pbc.donor.astype(str)
    S["pbc"] = pbc[pbc.compartment.isin(S["per"])].reset_index(drop=True)
    G = [None] * N_COMP
    t0 = time.time()
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, (i, g, R) in enumerate(pool.imap_unordered(companion_task, range(N_COMP))):
            G[i] = (g, R)
            if (k + 1) % 50 == 0:
                print(f"companions {k + 1}/{N_COMP} {time.time() - t0:.0f}s", flush=True)
    return G


# ---- power -----------------------------------------------------------------------------------------------------

def power(jobs):
    t0 = time.time()
    organs = setup()
    print(f"clock: {len(organs)} of {len(S['per'])} peripheral organs pass the precondition ({time.time() - t0:.0f}s)", flush=True)
    P, R, Z = output_sets(jobs)
    joint = (P & S["pass_c"][None, :]).sum(1)
    ok = joint >= MIN_JOINT
    print(f"output sets: passing organs median {np.median(P.sum(1)):.0f} (min {P.sum(1).min()}, max {P.sum(1).max()}); "
          f"joint with the clock: median {np.median(joint):.0f}; sets with >= {MIN_JOINT} joint organs: {int(ok.sum())}", flush=True)
    outs = out_locals(P, R, Z)
    G = companions(jobs)
    rows = []
    for i, (g, Rp) in enumerate(G):
        lp = standardise(local(Rp / S["Zc"]))
        m, n = rho_med(lp, outs)
        rows.append(dict(companion=i, rho_med=m, sets=n, genes=";".join(g)))
    C = pd.DataFrame(rows)
    C.to_csv(OUT / "companions.csv", index=False)
    q95 = float(C.rho_med.quantile(0.95))
    print(f"companions: rho_med median {C.rho_med.median():.4f}, q95 {q95:.4f}, q05 {C.rho_med.quantile(0.05):.4f} "
          f"({time.time() - t0:.0f}s)", flush=True)
    # injection: output residual = zhat_j * f * dl_c + own residuals permuted between animals within organ
    dlc = np.nan_to_num(local(S["Rc"] / S["Zc"]))
    rng = np.random.default_rng([SEED, 3])
    lc = standardise(local(S["Rc"] / S["Zc"]))
    levels = (0.0, 0.5, 1.0)
    res = {f: [] for f in levels}
    sel = np.flatnonzero(ok)
    for d in range(N_DRAW):
        Rp = {j: R[j].copy() for j in sel}
        for j in sel:
            for k in np.flatnonzero(P[j]):
                idx = np.flatnonzero(np.isfinite(R[j][:, k]))
                Rp[j][idx, k] = R[j][rng.permutation(idx), k]
        for f in levels:
            inj = []
            for j in sel:
                D = Rp[j] / Z[j]
                jt = P[j] & S["pass_c"]
                D = D + np.where(jt[None, :] & np.isfinite(D), f * dlc, 0)
                inj.append(standardise(local(D)))
            res[f].append(rho_med(lc, inj)[0])
        if (d + 1) % 20 == 0:
            print(f"injection {d + 1}/{N_DRAW} {time.time() - t0:.0f}s", flush=True)
    inj = pd.DataFrame([dict(inheritance=f, draws=N_DRAW, rho_med_median=float(np.median(v)), rho_med_q05=float(np.quantile(v, 0.05)),
                             power=float(np.mean(np.array(v) > q95))) for f, v in res.items()])
    inj.to_csv(OUT / "power_injection.csv", index=False)
    full = float(inj[inj.inheritance == 1.0].power.iloc[0])
    Pw = pd.DataFrame([dict(clock_organs=len(organs), peripheral_organs=len(S["per"]), sets=N_SETS, sets_ge15_joint=int(ok.sum()),
                            joint_median=float(np.median(joint)), companion_rho_med_median=float(C.rho_med.median()),
                            companion_rho_med_q95=q95, power_full=full,
                            power_half=float(inj[inj.inheritance == 0.5].power.iloc[0]),
                            power_zero=float(inj[inj.inheritance == 0.0].power.iloc[0]), decidable=bool(full >= 0.8),
                            runtime_s=round(time.time() - t0))])
    Pw.to_csv(OUT / "power.csv", index=False)
    print(inj.round(4).to_string(index=False))
    print(Pw.T.to_string(), flush=True)


# ---- score (not run before the gate is read) ---------------------------------------------------------------------

def score(jobs):
    organs = setup()
    P, R, Z = output_sets(jobs)
    ok = (P & S["pass_c"][None, :]).sum(1) >= MIN_JOINT
    outs = out_locals(P, R, Z)
    C = pd.read_csv(OUT / "companions.csv")
    q95 = float(C.rho_med.quantile(0.95))
    Dc = S["Rc"] / S["Zc"]
    lc = standardise(local(Dc))
    rows = []
    sel = np.flatnonzero(ok)
    rj = np.array([rho(lc, outs[j]) for j in sel])
    obs, n = float(np.nanmedian(rj)), len(sel)
    # permutation of the clock residuals between animals within organ (recipient's zhat)
    rng = np.random.default_rng([SEED, 4])
    null = []
    for _ in range(N_PERM):
        Rp = S["Rc"].copy()
        for k in np.flatnonzero(S["pass_c"]):
            idx = np.flatnonzero(np.isfinite(S["Rc"][:, k]))
            Rp[idx, k] = S["Rc"][rng.permutation(idx), k]
        null.append(rho_med(standardise(local(Rp / S["Zc"])), outs)[0])
    p = (1 + np.sum(np.array(null) >= obs)) / (1 + N_PERM)
    rows.append(dict(readout="H413.1 rho_med", value=obs, sets=n, companion_q95=q95, above_q95=bool(obs > q95), perm_p=p,
                     passes=bool(obs > q95 and p < 0.05)))
    for part in ("re", "im"):  # amplitude-only and phase-only parts
        v = [rho(lc, outs[j], part) for j in sel]
        rows.append(dict(readout=f"rho_med {part}", value=float(np.nanmedian(v)), sets=n))
    # shared part: animal-mean departures over passing organs, 12 animals
    mc = np.nanmean(Dc, axis=1)
    sh = [rho(mc[:, None], np.nanmean(R[j] / Z[j], axis=1)[:, None]) for j in sel]
    rows.append(dict(readout="shared part rho_med", value=float(np.nanmedian(sh)), sets=n))
    # by output phase bin
    sets, pool = pd.read_csv(M.OUT / "sets.csv"), pd.read_csv(M.OUT / "pool.csv").set_index("gene")
    binj = sets.groupby("set").gene.apply(lambda g: int(pool.loc[g, "bin"].mode().iloc[0]))
    for b in sorted(binj.unique()):
        js = [j for j in sel if binj[j] == b]
        if js:
            rows.append(dict(readout=f"rho_med phase bin {b}", value=float(np.nanmedian([rho(lc, outs[j]) for j in js])), sets=len(js)))
    # core loop: clock without the trio (same passing organs), and 200 random drops of 3 clock genes (gene-drop null)
    clock = list(WF.genes)
    def drop_read(drop):
        keep = [g for g in clock if g not in drop]
        ds = cs.Dataset(S["pb"][KEYS + keep], anchors=S["per"], waveform=S["wf_c"].subset(keep), name="a413 drop", min_donors=5)
        r, z = grids(ds, S["zt"], S["per"], S["animals"], organs)
        return rho_med(standardise(local(r / S["Zc"])), outs)[0]
    rows.append(dict(readout="core loop (trio dropped) rho_med", value=drop_read(list(TRIO))))
    rg = np.random.default_rng([SEED, 5])
    dn = [drop_read([clock[i] for i in rg.choice(len(clock), 3, replace=False)]) for _ in range(200)]
    rows.append(dict(readout="random 3-gene drops rho_med q05 / median", value=float(np.nanmedian(dn)), companion_q95=float(np.nanquantile(dn, 0.05))))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "tests.csv", index=False)
    print(T.round(4).to_string(index=False))


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    (power if a.stage == "power" else score)(min(a.jobs, 4))


if __name__ == "__main__":
    main()
