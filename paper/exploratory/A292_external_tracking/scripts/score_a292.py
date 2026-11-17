"""A292 (registered in `paper/docs/registered_experiments.md`): do weak clocks depart from the body's time toward the
external time of day?

GTEx fast deaths, the A104 rows (tissues with >= 40 donors, every tissue an anchor). Per periphery tissue k: clock-plane
states x of its samples (the fit's states), its leave-one-out host h, and u = exp(i 2 pi tod / 24) rotated into the
host frame (delta = arg sum h_full conj(u)). Complex least squares x = mu_stratum + a h + b u; tracking weight
w = Re(b / a) / sd(h). H292.1: Spearman rho(own amplitude A, w) < 0, one-sided p against the hour of death permuted
over donors (2000 permutations). b is obtained by partialling out [intercepts, h] (Frisch-Waugh), so all permutations
of one tissue are one matrix product.

    python score_a292.py --stage calibrate [--n-sim 500] [--n-perm 2000] [--jobs 4]   # V292, no observed statistic
    python score_a292.py --stage score [--n-sets 100] [--jobs 4]                      # H292.1, companions, descriptive
Outputs (`../output/`): calibration.csv; tissues.csv, tests.csv, sets.csv, descriptive.csv. Aggregates only (GTEx
time of death is controlled access: no per-donor value is written).
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync._core import host as H  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 292
BAND = 0.0685
META = ["donor", "compartment", "depth"]


# ---- design per tissue ------------------------------------------------------------------------------------------

def _cx(f):
    return f[:, 0] + 1j * f[:, 1]


def designs(ds: cs.Dataset, tod: pd.Series) -> dict:
    """Per compartment: states x, LOO host h (with posterior cov), donor ids, group dummies; and the full host."""
    d = ds.data
    comps = list(ds.anchors)
    hosts = H._hosts(d, comps, comps)
    full = hosts.pop("__full__")
    fd = np.array([str(r).split("#")[0] for r in full.resections])
    out = {"__full__": dict(h=_cx(full.f), donor=fd)}
    for k in comps:
        h = hosts[k]
        pos = pd.Series(np.arange(len(h.resections)), index=[str(r) for r in h.resections])
        i = np.flatnonzero((d.compartment == k) & np.isin(d.resection.astype(str), pos.index))
        donor = np.array([str(r).split("#")[0] for r in d.resection[i]])
        ok = tod.reindex(donor).notna().to_numpy()
        i, donor = i[ok], donor[ok]
        j = pos[d.resection[i].astype(str)].to_numpy()
        G = pd.get_dummies(pd.Series(d.group[i])).to_numpy(float)
        out[k] = dict(x=_cx(d.beta[i]), V=d.V[i], h=_cx(h.f[j]), hcov=h.f_cov[j], donor=donor, G=G)
    return out


def _proj(M):
    """Residual-maker I - M (M^H M)^-1 M^H."""
    Q, _ = np.linalg.qr(M)
    return np.eye(len(M)) - Q @ Q.conj().T


def weights(D: dict, comps, U: np.ndarray, x_over=None) -> np.ndarray:
    """w per compartment (rows) for each column of U (donor-indexed external-time vectors, already rotated into the
    host frame): w = Re(b / a) / sd(h). U: DataFrame-like dict donor -> row index via D['__index__']."""
    idx = D["__index__"]
    W = np.full((len(comps), U.shape[1]), np.nan)
    for r, k in enumerate(comps):
        e = D[k]
        x = e["x"] if x_over is None else x_over[k]
        Z = np.column_stack([e["G"], e["h"]])
        R = _proj(Z)
        u = U[idx[e["donor"]]]  # n x P
        ru, rx = R @ u, R @ x
        b = (ru.conj() * rx[:, None]).sum(0) / (np.abs(ru) ** 2).sum(0)
        coef, *_ = np.linalg.lstsq(Z, x[:, None] - u * b[None, :], rcond=None)  # a for each column
        a = coef[-1]
        hc = e["h"] - e["h"].mean()
        W[r] = np.real(b / a) / np.sqrt(np.mean(np.abs(hc) ** 2))
    return W


def external(D: dict, tod: pd.Series, perms: np.ndarray) -> np.ndarray:
    """u vectors for the real hours (column 0) and each permutation of the hours over donors, each rotated into the host
    frame by delta = arg sum_d h_full conj(u_d) (recomputed per column). Rows indexed by D['__index__']."""
    donors = D["__donors__"]
    t = tod.reindex(donors).to_numpy(float)
    T = np.column_stack([t] + [t[p] for p in perms]) if len(perms) else t[:, None]
    U = np.exp(1j * 2 * np.pi * T / 24)
    full = D["__full__"]
    fi = D["__index__"][full["donor"]]
    delta = np.angle((full["h"][:, None] * U[fi].conj()).sum(0))
    return U * np.exp(1j * delta)[None, :]


def index_donors(D, comps):
    donors = sorted(set(np.concatenate([D[k]["donor"] for k in comps] + [D["__full__"]["donor"]])))
    D["__donors__"] = np.array(donors)
    D["__index__"] = pd.Series(np.arange(len(donors)), index=donors)
    return D


def rho_lt(A, W):
    """Spearman rho(A, w) for each column of W (tissues with finite w)."""
    out = np.full(W.shape[1], np.nan)
    for c in range(W.shape[1]):
        ok = np.isfinite(W[:, c]) & np.isfinite(A)
        out[c] = spearmanr(A[ok], W[ok, c]).statistic
    return out


def permutations(n_donors, n_perm, seed):
    rng = np.random.default_rng([seed, 1])
    return np.array([rng.permutation(n_donors) for _ in range(n_perm)]) if n_perm else np.empty((0, n_donors), int)


def test(D, comps, A, tod, n_perm, seed):
    perms = permutations(len(D["__donors__"]), n_perm, seed)
    U = external(D, tod, perms)
    W = weights(D, comps, U)
    r = rho_lt(A, W)
    p = (1 + np.sum(r[1:] <= r[0])) / (1 + n_perm)
    return float(r[0]), float(p), W[:, 0]


# ---- calibration (V292) -----------------------------------------------------------------------------------------

def simulate_states(D, comps, tod, fit, rng, w0=0.0, A=None):
    """S2 (w0 = 0) or S3 states at the real design: truth = full host; tissue states a h_true + own noise + measurement
    noise (+ b u_real for S3); LOO host = truth + posterior noise. Returns (x_over, D_sim) with simulated hosts."""
    full = D["__full__"]
    h_true = pd.Series(full["h"], index=full["donor"])
    U0 = external(D, tod, np.empty((0, 0), int))[:, 0]
    Ds = {kk: vv for kk, vv in D.items()}
    xs = {}
    amin = np.nanmin(A)
    for r, k in enumerate(comps):
        e = D[k]
        ht = h_true.reindex(e["donor"]).to_numpy()
        ok = np.isfinite(ht)
        a, mu, psi = fit[k]
        n = len(e["donor"])
        noise = rng.standard_normal((n, 2)) * np.sqrt(psi)
        L = np.linalg.cholesky(e["V"] + 1e-12 * np.eye(2))
        noise += np.einsum("nab,nb->na", L, rng.standard_normal((n, 2)))
        x = mu + a * np.where(ok, ht, 0) + noise[:, 0] + 1j * noise[:, 1]
        if w0:
            hc = e["h"] - e["h"].mean()
            b = w0 * a * np.sqrt(np.mean(np.abs(hc) ** 2)) * amin / A[r]
            x = x + b * U0[D["__index__"][e["donor"]]]
        Lh = np.linalg.cholesky(e["hcov"] + 1e-12 * np.eye(2))
        hn = np.einsum("nab,nb->na", Lh, rng.standard_normal((n, 2)))
        Ds[k] = dict(e, h=np.where(ok, ht, e["h"]) + hn[:, 0] + 1j * hn[:, 1])
        xs[k] = x
    return xs, Ds


def tissue_fit(D, comps, ds_table):
    """a (on the LOO host), mean mu and own variance psi (A104 table, per dimension) per tissue."""
    out = {}
    for k in comps:
        e = D[k]
        Z = np.column_stack([np.ones(len(e["x"])), e["h"]])
        coef, *_ = np.linalg.lstsq(Z, e["x"], rcond=None)
        out[k] = (coef[1], coef[0], float(ds_table.loc[k, "psi"]))
    return out


def one_sim(D, comps, tod, fit, A, w0, seed, n_perm):
    limit_threads(1)
    rng = np.random.default_rng([SEED, 7, int(w0 * 100), seed])
    xs, Ds = simulate_states(D, comps, tod, fit, rng, w0, A)
    perms = permutations(len(D["__donors__"]), n_perm, SEED + 1000 + seed)
    U = external(Ds, tod, perms)
    W = weights(Ds, comps, U, x_over=xs)
    r = rho_lt(A, W)
    return float((1 + np.sum(r[1:] <= r[0])) / (1 + n_perm))


# ---- main ---------------------------------------------------------------------------------------------------------

def block_dataset(pb, block):
    comps = sorted(c for c in pb.compartment.unique() if c.startswith("Brain") == (block == "brain"))
    return cs.Dataset(pb[pb.compartment.isin(comps)], anchors=comps, name=f"gtex_{block}")


def prepared(pb, tod, block):
    ds = block_dataset(pb, block)
    pt = A104.per_tissue(ds)
    comps = [k for k in ds.anchors if k in pt.index]
    D = index_donors(designs(ds, tod), comps)
    return ds, pt, comps, D


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("calibrate", "score"), required=True)
    ap.add_argument("--n-sim", type=int, default=500)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--n-sets", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    warnings.filterwarnings("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    pb, donors = gtex.load(min_donors=A104.MIN_DONORS)
    limit_threads(1)
    tod = donors.tod
    ds, pt, comps, D = prepared(pb, tod, "periphery")
    A = pt.loc[comps, "own_amplitude"].to_numpy(float)
    rep = A104.rho(pt)
    print(f"periphery: {len(comps)} tissues, {len(D['__donors__'])} donors; A104 rho {rep:.3f} {time.time() - t0:.0f}s", flush=True)
    if round(rep, 3) != 0.734:
        raise SystemExit("A104 not reproduced")

    if a.stage == "calibrate":
        fit = tissue_fit(D, comps, pt)
        n_perm = 500  # inside each simulation (stated in the RP entry); the scored test uses --n-perm
        rows = []
        for w0 in (0.0, 0.1, 0.25):
            ps = Parallel(n_jobs=a.jobs)(delayed(one_sim)(D, comps, tod, fit, A, w0, s, n_perm) for s in range(a.n_sim))
            ps = np.array(ps)
            rows.append(dict(scenario="S2 (no tracking)" if w0 == 0 else f"S3 w0={w0}", w0=w0, n_sim=len(ps),
                             n_perm_inside=n_perm, rejection=float(np.mean(ps < 0.05))))
            print(rows[-1], f"{time.time() - t0:.0f}s", flush=True)
        cal = pd.DataFrame(rows)
        cal["V292_calibrated"] = cal.rejection.where(cal.w0 == 0).le(BAND).where(cal.w0 == 0)
        cal.to_csv(OUT / "calibration.csv", index=False)
        print(cal.to_string(), flush=True)
        return

    cal = pd.read_csv(OUT / "calibration.csv")
    if not bool(cal.loc[cal.w0 == 0, "rejection"].iloc[0] <= BAND):
        raise SystemExit("V292 failed: H292.1 is not read (registered)")
    r0, p, w = test(D, comps, A, tod, a.n_perm, SEED)
    power25 = float(cal.loc[cal.w0 == 0.25, "rejection"].iloc[0])
    tis = pd.DataFrame(dict(block="periphery", tissue=comps, own_amplitude=A, coupling=pt.loc[comps, "coupling"].to_numpy(),
                            w=w, n=[len(D[k]["donor"]) for k in comps]))
    rows = [dict(block="periphery", analysis="H292.1", tissues=len(comps), rho=r0, p_one_sided=p,
                 mean_w=float(np.mean(w)), median_w=float(np.median(w)), power_w0_025=power25,
                 passes=bool(r0 < 0 and p < 0.05))]
    print(rows[-1], f"{time.time() - t0:.0f}s", flush=True)

    # descriptive: brain block, ventilated donors (periphery and brain)
    desc = []
    for label, hardy in (("fast", gtex.FAST), ("ventilated", (0,))):
        pbx, dx = (pb, donors) if label == "fast" else gtex.load(min_donors=A104.MIN_DONORS, hardy=hardy)
        for block in ("periphery", "brain"):
            if label == "fast" and block == "periphery":
                continue
            try:
                _, ptx, cx, Dx = prepared(pbx, dx.tod, block)
                Ax = ptx.loc[cx, "own_amplitude"].to_numpy(float)
                rx, px, wx = test(Dx, cx, Ax, dx.tod, a.n_perm, SEED + 1)
                desc.append(dict(deaths=label, block=block, tissues=len(cx), rho=rx, p_one_sided=px,
                                 median_w=float(np.median(wx)), mean_w=float(np.mean(wx))))
                tis = pd.concat([tis, pd.DataFrame(dict(block=f"{block}|{label}", tissue=cx, own_amplitude=Ax,
                                                        coupling=ptx.loc[cx, "coupling"].to_numpy(), w=wx,
                                                        n=[len(Dx[k]["donor"]) for k in cx]))])
            except Exception as e:  # noqa: BLE001
                desc.append(dict(deaths=label, block=block, error=f"{type(e).__name__}: {e}"[:200]))
            print(desc[-1], f"{time.time() - t0:.0f}s", flush=True)
    pd.DataFrame(desc).to_csv(OUT / "descriptive.csv", index=False)
    tis.to_csv(OUT / "tissues.csv", index=False)

    # companions: the statistic (no permutation) for matched and responsive pseudo-clock sets
    sets = companions(pb_genes(), tod, a.n_sets, a.jobs)
    sets.to_csv(OUT / "sets.csv", index=False)
    for kind in ("matched", "responsive"):
        s = sets[(sets.kind == kind)].rho.dropna()
        rows[0][f"{kind}_median"] = float(s.median())
        rows[0][f"{kind}_pi_below"] = float(np.mean(s > r0))  # share of sets with a less negative rho
        rows[0][f"{kind}_n"] = int(len(s))
    pd.DataFrame(rows).to_csv(OUT / "tests.csv", index=False)
    print(pd.DataFrame(rows).T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


def pb_genes():
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(a.var_names)
    a.file.close()
    return gtex.load(min_donors=A104.MIN_DONORS, genes=genes)


S = {}


def _set_task(i):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    kind, w = S["sets"][i]
    pb, tod = S["pb"], S["tod"]
    comps = sorted(c for c in pb.compartment.unique() if not c.startswith("Brain"))
    try:
        ds = cs.Dataset(pb.loc[pb.compartment.isin(comps), META + list(w.genes)], anchors=comps, waveform=w, name=f"s{i}")
        pt = A104.per_tissue(ds)
        cx = [k for k in ds.anchors if k in pt.index]
        D = index_donors(designs(ds, tod), cx)
        U = external(D, tod, np.empty((0, 0), int))
        W = weights(D, cx, U)
        r = rho_lt(pt.loc[cx, "own_amplitude"].to_numpy(float), W)[0]
        return dict(set=i, kind=kind, rho=float(r))
    except Exception as e:  # noqa: BLE001
        return dict(set=i, kind=kind, error=f"{type(e).__name__}: {e}"[:200])


def companions(loaded, tod, n_sets, jobs):
    from multiprocessing import Pool
    from clocksync.paired import pseudo_waveforms
    pb, _ = loaded
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    del x
    wf = cs.gtex()
    sets = [("matched", w) for w in pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, 1]), wf, kind="matched")]
    sets += [("responsive", w) for w in pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, 2]), wf, kind="responsive")]
    keep = sorted(set(g for _, w in sets for g in w.genes))
    S.update(pb=pb[META + keep].copy(), tod=tod, sets=sets)
    with Pool(jobs) as pool:
        return pd.DataFrame(pool.map(_set_task, range(len(sets))))


if __name__ == "__main__":
    main()
