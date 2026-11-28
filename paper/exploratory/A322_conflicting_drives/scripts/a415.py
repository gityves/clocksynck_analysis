"""A415 (registered in `paper/docs/registered_experiments.md`): when conflicting drives cost a mouse organ clock
amplitude, does the organ's own clock noise psi stay (uptake) or rise (animal-specific split between the drives)?

Seven timed units: Manella 2021 Alb-Cre controls (ad libitum vs daytime feeding; Liver, WAT, Lung) and Xin 2021 (NRF vs
DRF; liver, heart, kidney, VAT). Per unit: clock states of A297 (`A297.liver_states`: one projector per tissue over all
its libraries, QC'd clock genes). Per arm: timed amplitude |a| (A297 OLS on ZT harmonics) and own noise psi = ML
residual variance of the state around the arm's timed fit beyond the Poisson measurement covariance of each library
(`clocksync._core.factor._known_f`, the estimator of `Dataset.timed_loadings`, applied to the states with their full
2 x 2 covariances). dA = log |a|^2 control - log |a|^2 conflict, dpsi = log psi control - log psi conflict (control:
AlbCre_AL / NRF; conflict: AlbCre_DF / DRF). Companion excess dpsi* = dpsi - median over 100 expression-matched
non-clock gene sets of the dataset (A263 rule: `pseudo_waveforms`, mean log fraction over all libraries, genes with mean
count >= 1; read at the clock positions, same libraries, same arms). kappa = median_u dpsi*_u / median_u dA_u.

  --stage power  simulation G415.1 (state-level model at each arm's design and measurement covariances; amplitudes from
                 the data, control psi from the data, conflict psi by truth A / B / C); companions of the units.
                 Prints counts and gate values only, never a psi or dpsi of the real data.
  --stage score  (written, NOT run before the power note) kappa with a 90 % interval (2,000 library bootstraps within
                 ZT and arm; the companion median is recomputed in a 50-set subsample per draw), per-unit table.

    paper/common/run.sh paper/exploratory/A322_conflicting_drives/scripts/a415.py --stage power|score [--check]
Outputs (`../output/a415/`): companions.csv, companions_summary.csv, power_sim.csv, power_summary.csv; score_*.csv.
"""
import argparse
import importlib.util
import re
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core import factor  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import FLOOR, projector, pseudo_waveforms, states  # noqa: E402
from common import data, mouse_ids  # noqa: E402
from common.paths import DATA  # noqa: E402

spec = importlib.util.spec_from_file_location("a415_a297", PAPER / "exploratory/A297_feeding_model/scripts/score_a297.py")
A297 = importlib.util.module_from_spec(spec)
sys.modules["a415_a297"] = A297
spec.loader.exec_module(A297)
QC = A297.C  # A111 `_common`: gene_qc, waveform_without

OUT = HERE.parents[1] / "output" / "a415"
W = 2 * np.pi / 24
WF = cs.gtex()
SEED = 415
N_SETS, N_SUB, N_BOOT, N_REP = 100, 50, 2000, 100
UNITS = (("manella2021", "Liver", "AlbCre_AL", "AlbCre_DF"), ("manella2021", "WAT", "AlbCre_AL", "AlbCre_DF"),
         ("manella2021", "Lung", "AlbCre_AL", "AlbCre_DF"), ("xin2021", "liver", "NRF", "DRF"),
         ("xin2021", "heart", "NRF", "DRF"), ("xin2021", "kidney", "NRF", "DRF"), ("xin2021", "VAT", "NRF", "DRF"))
XIN_N_READS = 30e6


# ---- states and estimators ---------------------------------------------------------------------------------------

def states_cov(x: pd.DataFrame, wf, depth):
    """Clock states beta (n x 2, as `A297.liver_states`: one projector over all rows of x) and the full Poisson
    measurement covariance V (n x 2 x 2) of each library."""
    L = projector(x[list(wf.genes)], wf, counts=True)
    beta, _ = states(x, L, wf, depth=depth)
    var = 1.0 / np.maximum(x[list(wf.genes)].to_numpy(float), FLOOR)
    return beta, np.einsum("ag,ng,bg->nab", L, var, L)


def amp_ols(z, zt) -> complex:
    """A297 `amplitude`: complex first-harmonic coefficient of the complex state z at ZT zt (h)."""
    t = np.asarray(zt, float) * W
    D = np.c_[np.ones(len(t)), np.exp(1j * t), np.exp(-1j * t)]
    return complex(np.linalg.lstsq(D, np.asarray(z), rcond=None)[0][1])


def timed_psi(beta, V, zt):
    """Timed fit of one arm (one intercept): (z ML amplitude, psi own noise beyond measurement noise)."""
    t = np.asarray(zt, float) * W
    f = np.column_stack([np.cos(t), np.sin(t)])
    _, _, z, psi, _ = factor._known_f(beta, V, f, np.zeros(len(t), int))
    return complex(z[0]), float(psi[0])


def arm_stats(beta, V, zt):
    """(|a_ols|^2, psi) of one arm."""
    a = amp_ols(beta[:, 0] + 1j * beta[:, 1], zt)
    return abs(a) ** 2, timed_psi(beta, V, zt)[1]


def lratio(a, b):
    return float(np.log(a / b)) if np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0 else np.nan


def kappa(dpsi, dA, comp=0.0):
    num, den = np.nanmedian(np.asarray(dpsi, float) - comp), np.nanmedian(dA)
    return float(num / den) if np.isfinite(num) and np.isfinite(den) and den != 0 else np.nan


# ---- units ------------------------------------------------------------------------------------------------------

def load_units():
    """Per unit: dict(x counts table of the tissue (all arms), beta, V, zt, arm, wf, idx0, idx1). The real states
    are checked against `A297.liver_states` (same projector, same genes)."""
    out, qcs = [], {}
    for name in ("manella2021", "xin2021"):
        pb, donors = data.load(name)
        qc = QC.gene_qc(pb, name)
        qcs[name] = qc
        wf = QC.waveform_without(sorted(qc[~qc.kept].gene.unique()))
        for nm, tissue, ctl, con in (u for u in UNITS if u[0] == name):
            x = pb[pb.compartment == tissue].set_index("donor")
            beta, V = states_cov(x, wf, x.depth.to_numpy(float))
            d = donors.reindex(x.index)
            ref, _ = A297.liver_states(name, tissue, "arm")
            assert np.allclose(ref.z.to_numpy(), beta[:, 0] + 1j * beta[:, 1]) and (ref.index == x.index).all()
            arm, zt = d.arm.astype(str).to_numpy(), d.zt.to_numpy(float)
            out.append(dict(dataset=name, tissue=tissue, ctl=ctl, con=con, beta=beta, V=V, zt=zt, arm=arm, wf=wf,
                            depth=x.depth.to_numpy(float), order=x.index.astype(str).to_numpy(), n=pd.Series(arm).value_counts().to_dict()))
    return out, qcs


def arm_rows(u, which):
    return np.flatnonzero(u["arm"] == u[which])


# ---- genome-wide tables for the companions ------------------------------------------------------------------------

def genome_manella():
    """{tissue: (counts libraries x genes (upper-case symbols), depth)} of all four arms (as the projector of
    `A297.liver_states`); depth as the builder (sum over the authors' retained genes)."""
    from common.datasets.manella2021 import build as M
    out = {}
    for t in M.TISSUES:
        x = pd.concat([M.read_tissue(t, g) for g in M.GENOTYPES], axis=1)
        out[t] = _tidy(x)
    return out


def genome_xin():
    from common.datasets.xin2021 import build as X
    out = {}
    for t, (f, how) in X.FILES.items():
        x = pd.read_csv(X.SRC / f, sep="\t", dtype=str, skiprows=[1])
        ids = x["geneID"] if "geneID" in x else x.iloc[:, -1]
        cols = [c for c in x.columns if re.fullmatch(r"[DN]\d+_\d+", c)]
        tpm = x[cols].astype(float)
        tpm.index = ids.to_numpy()
        tpm.columns = [f"{'DRF' if c[0] == 'D' else 'NRF'}_ZT{int(re.search(r'[DN](\d+)_', c)[1])}_m{c.split('_')[1]}"
                       for c in cols]
        sym = mouse_ids.to_symbols(tpm.index, how)
        tpm = tpm[sym.notna().to_numpy()]
        tpm.index = sym[sym.notna()].to_numpy()
        out[t] = _tidy((tpm * XIN_N_READS / 1e6).round(), depth=XIN_N_READS)
    return out


def _tidy(x: pd.DataFrame, depth=None):
    """genes x libraries -> (libraries x canonical symbols, depth); duplicated symbols are summed."""
    depth = x.sum() if depth is None else pd.Series(float(depth), index=x.columns)
    x = x.copy()
    x.index = [canonical(g) for g in x.index]
    x = x.groupby(level=0).sum()
    return x.T, depth


def companions(units, genome, rng_tag, t0):
    """dpsi and dA of 100 expression-matched sets per unit (rows: dataset, tissue, set). Returns (table, summary)."""
    rows, summ = [], []
    for name in ("manella2021", "xin2021"):
        g = genome[name]
        allx = pd.concat([x for x, _ in g.values()])
        alld = pd.concat([d for _, d in g.values()]).reindex(allx.index)
        lm = pd.Series(np.log((allx.to_numpy(np.float32) / alld.to_numpy(float)[:, None]).mean(0) + 1e-12), index=allx.columns)
        lm = lm[allx.mean(0).to_numpy() >= 1]
        miss = [c for c in WF.genes if c not in lm.index]
        if miss:  # clock genes the authors' filter removed: neutral expression for the bin, the real QC drops them anyway
            lm = pd.concat([lm, pd.Series(float(lm.median()), index=miss)])
        sets = pseudo_waveforms(lm, N_SETS, np.random.default_rng([SEED, rng_tag[name]]), WF, "matched")
        print(f"companions {name}: {len(lm)} matching genes, {len(sets)} sets ({time.time() - t0:.0f}s)", flush=True)
        for u in (u for u in units if u["dataset"] == name):
            x, depth = g[u["tissue"]]
            x, depth = x.reindex(u["order"]), depth.reindex(u["order"])  # libraries in the real pseudobulk's order
            assert not x.isna().all(axis=1).any()
            i0, i1 = arm_rows(u, "ctl"), arm_rows(u, "con")
            # QC (A111): a set position is dropped in a tissue of the dataset if its median count < 5
            fails = {}
            for tk, (xt, _) in g.items():
                mt = xt.median()
                for w in sets:
                    fails.setdefault(w.name, set()).update(c for c in w.genes if mt.get(c, 0) < QC.EXPR_MIN)
            for w in sets:
                try:
                    wk = w.drop(sorted(fails[w.name])) if fails[w.name] else w
                    xs = x[list(wk.genes)]
                    beta, V = states_cov(xs, wk, depth.to_numpy(float))
                    a0, p0 = arm_stats(beta[i0], V[i0], u["zt"][i0])
                    a1, p1 = arm_stats(beta[i1], V[i1], u["zt"][i1])
                    rows.append(dict(dataset=name, tissue=u["tissue"], set=w.name, n_genes=len(wk.genes),
                                     dpsi=lratio(p0, p1), dA=lratio(a0, a1), error=""))
                except Exception as e:
                    rows.append(dict(dataset=name, tissue=u["tissue"], set=w.name, n_genes=0, dpsi=np.nan, dA=np.nan,
                                     error=f"{type(e).__name__}: {e}"[:120]))
            print(f"companions {name} {u['tissue']} done ({time.time() - t0:.0f}s)", flush=True)
    c = pd.DataFrame(rows)
    s = c.groupby(["dataset", "tissue"]).agg(sets=("set", "nunique"), finite=("dpsi", lambda v: int(np.isfinite(v).sum())),
                                              failed=("error", lambda v: int((v != "").sum())),
                                              min_genes=("n_genes", lambda v: int(v[v > 0].min()) if (v > 0).any() else 0)).reset_index()
    return c, s


def comp_matrix(c):
    """sets x units matrix of companion dpsi, columns in UNITS order (tissue keys)."""
    c = c.assign(unit=c.dataset + ":" + c.tissue)
    m = c.pivot(index="set", columns="unit", values="dpsi")
    return m[[f"{d}:{t}" for d, t, _, _ in UNITS]]


# ---- G415.1 simulation -----------------------------------------------------------------------------------------

def sim_arm(u, rows, z, psi, rng):
    """State-level model beta = M(z)(cos, sin) + e, e ~ N(0, psi I + V_i) at the arm's design (ZT, V_i fixed)."""
    t = u["zt"][rows] * W
    s = z * np.exp(1j * t)
    mean = np.c_[s.real, s.imag]
    C = u["V"][rows] + max(psi, 1e-12) * np.eye(2)
    e = np.einsum("nab,nb->na", np.linalg.cholesky(C), rng.standard_normal((len(t), 2)))
    return mean + e


def power(a):
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    sfx = "_check" if a.check else ""
    units, _ = load_units()
    print("libraries per arm (control / conflict): " + "; ".join(
        f"{u['dataset']}:{u['tissue']} {u['n'][u['ctl']]}/{u['n'][u['con']]}" for u in units), flush=True)
    # (1) companions, if genome-wide counts exist
    genome = {}
    try:
        genome = {"manella2021": genome_manella(), "xin2021": genome_xin()}
    except Exception as e:
        print(f"no genome-wide counts: {type(e).__name__}: {e}", flush=True)
    have = bool(genome)
    if have:
        if a.check:
            global N_SETS
            N_SETS = 8
        c, s = companions(units, genome, {"manella2021": 1, "xin2021": 2}, t0)
        c.to_csv(OUT / f"companions{sfx}.csv", index=False)
        s.to_csv(OUT / f"companions_summary{sfx}.csv", index=False)
        cm = comp_matrix(c).to_numpy(float)
        print(s.to_string(index=False), flush=True)
    else:
        cm = np.zeros((1, len(units)))
        pd.DataFrame([dict(genome_wide="absent")]).to_csv(OUT / "companions_summary.csv", index=False)
    full = np.nanmedian(cm, axis=0)
    # (2) simulation inputs: amplitudes of both arms (ML), control psi; never written
    inp = []
    for u in units:
        i0, i1 = arm_rows(u, "ctl"), arm_rows(u, "con")
        z0, p0 = timed_psi(u["beta"][i0], u["V"][i0], u["zt"][i0])
        z1, p1 = timed_psi(u["beta"][i1], u["V"][i1], u["zt"][i1])
        inp.append((i0, i1, z0, p0, z1))
    n_rep = 5 if a.check else N_REP
    rows = []
    for truth in "ABC":
        rng = np.random.default_rng([SEED, ord(truth)])
        for r in range(n_rep):
            dp, dA = [], []
            for u, (i0, i1, z0, p0, z1) in zip(units, inp):
                psi1 = {"A": p0, "B": p0 * abs(z1) ** 2 / abs(z0) ** 2, "C": 2 * p0}[truth]
                b0, b1 = sim_arm(u, i0, z0, p0, rng), sim_arm(u, i1, z1, psi1, rng)
                s0, q0 = arm_stats(b0, u["V"][i0], u["zt"][i0])
                s1, q1 = arm_stats(b1, u["V"][i1], u["zt"][i1])
                dp.append(lratio(q0, q1)), dA.append(lratio(s0, s1))
            sub = rng.choice(len(cm), min(N_SUB, len(cm)), replace=False)  # companion median subsampling noise only
            noise = np.nanmedian(cm[sub], axis=0) - full
            rows.append(dict(truth=truth, rep=r, kappa=kappa(dp, dA, noise), median_dA=float(np.nanmedian(dA))))
        print(f"truth {truth} done ({time.time() - t0:.0f}s)", flush=True)
    sim = pd.DataFrame(rows)
    sim.to_csv(OUT / f"power_sim{sfx}.csv", index=False)
    k = {t: sim[sim.truth == t].kappa for t in "ABC"}
    sd_a = float(k["A"].std(ddof=1))
    share_a = float(((k["A"].abs() + 1.645 * sd_a) < 0.5).mean())
    summ = dict(units=len(units), companions=have, n_sets=N_SETS if have else 0, n_rep=n_rep,
                kappa_A_median=float(k["A"].median()), kappa_A_sd=sd_a, G415_1_A_share_inside=share_a,
                G415_1_A_pass=bool(share_a >= 0.8), kappa_B_median=float(k["B"].median()),
                G415_1_B_pass=bool(k["B"].median() >= 0.7), kappa_C_median=float(k["C"].median()),
                G415_1_C_pass=bool(k["C"].median() <= -0.3), n_nan_kappa=int(sim.kappa.isna().sum()))
    summ["G415_1_decidable"] = bool(summ["G415_1_A_pass"] and summ["G415_1_B_pass"] and summ["G415_1_C_pass"])
    pd.Series(summ).to_csv(OUT / f"power_summary{sfx}.csv", header=["value"])
    print(pd.Series(summ).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


# ---- score -----------------------------------------------------------------------------------------------------

def boot_unit(u, rng, n):
    """n library bootstraps within ZT and arm: arrays (dpsi, dA)."""
    groups = {w: [np.flatnonzero((u["arm"] == u[w]) & (u["zt"] == z)) for z in np.unique(u["zt"][u["arm"] == u[w]])]
              for w in ("ctl", "con")}
    dp, dA = np.full(n, np.nan), np.full(n, np.nan)
    for j in range(n):
        st = {}
        for w, gs in groups.items():
            i = np.concatenate([rng.choice(g, len(g)) for g in gs])
            st[w] = arm_stats(u["beta"][i], u["V"][i], u["zt"][i])
        dA[j], dp[j] = lratio(st["ctl"][0], st["con"][0]), lratio(st["ctl"][1], st["con"][1])
    return dp, dA


def score(a):
    limit_threads(4)
    t0 = time.time()
    ps = pd.read_csv(OUT / "power_summary.csv", index_col=0).value
    decidable = str(ps["G415_1_decidable"]) == "True"
    c = pd.read_csv(OUT / "companions.csv") if (OUT / "companions.csv").exists() else None
    units, _ = load_units()
    cm = comp_matrix(c).to_numpy(float) if c is not None else np.zeros((1, len(units)))
    cmed = np.nanmedian(cm, axis=0)
    rng = np.random.default_rng([SEED, 99])
    dp, dA, tab = [], [], []
    BD, BA = [], []
    for k, u in enumerate(units):
        i0, i1 = arm_rows(u, "ctl"), arm_rows(u, "con")
        s0, q0 = arm_stats(u["beta"][i0], u["V"][i0], u["zt"][i0])
        s1, q1 = arm_stats(u["beta"][i1], u["V"][i1], u["zt"][i1])
        dp.append(lratio(q0, q1)), dA.append(lratio(s0, s1))
        bd, ba = boot_unit(u, rng, 200 if a.check else N_BOOT)
        BD.append(bd), BA.append(ba)
        tab.append(dict(dataset=u["dataset"], tissue=u["tissue"], n_ctl=len(i0), n_con=len(i1), dA=dA[-1], dpsi=dp[-1],
                        comp_dpsi=cmed[k], dpsi_star=dp[-1] - cmed[k]))
        print(f"{u['tissue']} ({time.time() - t0:.0f}s)", flush=True)
    BD, BA = np.array(BD), np.array(BA)
    t = pd.DataFrame(tab)
    sel = {"all": np.ones(len(units), bool), **{d: np.array([u["dataset"] == d for u in units]) for d in ("manella2021", "xin2021")}}
    rows = []
    for lab, m in sel.items():
        k_obs = kappa(np.array(dp)[m], np.array(dA)[m], cmed[m])
        ks = []
        for j in range(BD.shape[1]):
            sub = rng.choice(len(cm), min(N_SUB, len(cm)), replace=False)
            ks.append(kappa(BD[m, j], BA[m, j], np.nanmedian(cm[sub][:, m], axis=0)))
        ks = np.array(ks)
        lo, hi = np.nanquantile(ks, [0.05, 0.95])
        verdict = ("pass" if -0.5 < lo and hi < 0.5 else "fail_below" if (hi < -0.5 or k_obs < -0.5) else
                   "fail_above" if (lo > 0.5 or k_obs > 0.5) else "inconclusive")
        rows.append(dict(subset=lab, units=int(m.sum()), kappa=k_obs, lo90=lo, hi90=hi, n_boot_ok=int(np.isfinite(ks).sum()),
                         H415_1=verdict if lab == "all" else "descriptive",
                         status="registered" if decidable and lab == "all" else "descriptive (G415.1 not met)" if lab == "all" else ""))
    out = pd.DataFrame(rows)
    t.to_csv(OUT / "score_units.csv", index=False)
    out.to_csv(OUT / "score_kappa.csv", index=False)
    pd.set_option("display.width", 200)
    print(t.round(3).to_string(index=False), out.round(3).to_string(index=False), sep="\n")


def main():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    (power if a.stage == "power" else score)(a)


if __name__ == "__main__":
    main()
