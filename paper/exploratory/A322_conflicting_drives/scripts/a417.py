"""A417 (registered in `paper/docs/registered_experiments.md`): does the fixed own noise psi mark perturbations of the
input (psi unchanged, kappa about 0) and its loss perturbations of the oscillator (psi follows the gain, kappa about
0.5-1)? kappa_arm = median_u dpsi*_u / median_u dA_u (A415 statistic); H417.1: kappa_O - kappa_I > 0.

Units (perturbed vs control, same laboratory, timed; psi beyond Poisson noise, `factor._known_f`, one fit per arm, on the
A297 clock states = one projector per tissue over all its libraries; companions A263 rule, 100 sets per dataset):
  O  Manella 2021 Liver, BLKO vs Alb-Cre, AL and DF (2 units)
  I  Greenwell 2019 liver arrhythmic vs ad libitum; GSE296764 liver and heart VIPKO vs WT (AL, TRF; CT0 / CT12:
     a = (z_CT0 - z_CT12) / 2, psi from the within-time-point scatter); Manella WAT and Lung BLKO vs Alb-Cre (AL, DF)
Control arm precondition (A395): timed |a| above the q95 of 1,000 time-label permutations of the control arm.
Implementation choices beyond the registration: see the docstrings of `amp` (2-level design), `boot_draw` (libraries
resampled within (arm, time point) by donor, the same donors in all tissues of a dataset) and `simulate` (replicate
interval from `N_BOOT_SIM` bootstraps; companions enter the draws only as the subsampling noise of the median).

  --stage power  QC, precondition, companions (output/a417/companions.csv), simulation G417.1 (truths S and N, 100
                 replicates). Prints counts and gate values only, never a psi, dpsi or perturbed-arm psi of the real data.
  --stage score  (written, NOT run before the power note) kappa_O, kappa_I, kappa_O - kappa_I with 90 % intervals
                 (2,000 bootstraps), H417.1-H417.3, per-unit / per-feeding / without-VIP-KO descriptives.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A322_conflicting_drives/scripts/a417.py --stage power|score [--check] [--keep-all]
Outputs (`../output/a417/`): power_units.csv, companions.csv, companions_summary.csv, power_sim.csv, power_summary.csv,
descriptive_loadability.csv; score_*.csv.
"""
import argparse
import importlib.util
import multiprocessing as mp
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync._core import factor  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import data, mouse_ids  # noqa: E402
from common.paths import DATA  # noqa: E402

spec = importlib.util.spec_from_file_location("a417_a415", HERE.parent / "a415.py")
A415 = importlib.util.module_from_spec(spec)
sys.modules["a417_a415"] = A415
spec.loader.exec_module(A415)
A297, QC, W, WF = A415.A297, A415.QC, A415.W, A415.WF
states_cov, lratio, sim_arm = A415.states_cov, A415.lratio, A415.sim_arm

OUT = HERE.parents[1] / "output" / "a417"
SEED = 417
N_SETS, N_SUB, N_BOOT, N_REP, N_BOOT_SIM, N_PERM = 100, 50, 2000, 100, 300, 1000
WORKERS = 4

# dataset -> (tissues, arm column, time column, QC compartments)
DATASETS = {"manella2021": (("Liver", "WAT", "Lung"), "arm", "zt", None),
            "greenwell2019": (("liver",), "arm", "zt", None),
            "vipko_gse296764": (("liver", "heart"), "arm", "ct", ("liver", "heart"))}
# (dataset, tissue, group, control arm, perturbed arm, feeding)
UNITS = (("manella2021", "Liver", "O", "AlbCre_AL", "BLKO_AL", "AL"), ("manella2021", "Liver", "O", "AlbCre_DF", "BLKO_DF", "DF"),
         ("greenwell2019", "liver", "I", "ad_libitum", "arrhythmic", "AL"),
         ("vipko_gse296764", "liver", "I", "WT_AL", "VIPKO_AL", "AL"), ("vipko_gse296764", "liver", "I", "WT_TRF", "VIPKO_TRF", "TRF"),
         ("vipko_gse296764", "heart", "I", "WT_AL", "VIPKO_AL", "AL"), ("vipko_gse296764", "heart", "I", "WT_TRF", "VIPKO_TRF", "TRF"),
         ("manella2021", "WAT", "I", "AlbCre_AL", "BLKO_AL", "AL"), ("manella2021", "WAT", "I", "AlbCre_DF", "BLKO_DF", "DF"),
         ("manella2021", "Lung", "I", "AlbCre_AL", "BLKO_AL", "AL"), ("manella2021", "Lung", "I", "AlbCre_DF", "BLKO_DF", "DF"))


def ulabel(u):
    return f"{u['dataset']}:{u['tissue']}:{u['ctl']}>{u['con']}"


# ---- estimators ----------------------------------------------------------------------------------------------------

def amp(z, zt) -> complex:
    """Complex first-harmonic amplitude of the state z at times zt (h): A297 OLS (mu + A e^{iwt} + B e^{-iwt}); with two
    time levels (VIP-KO CT0 / CT12, where e^{iwt} = +-1 and e^{-iwt} is collinear) the design is mu + A e^{iwt}, so
    A = (mean z_CT0 - mean z_CT12) / 2 (the CT0 phase convention)."""
    t = np.asarray(zt, float) * W
    cols = [np.ones(len(t)), np.exp(1j * t)] + ([np.exp(-1j * t)] if len(np.unique(np.round(t % (2 * np.pi), 6))) > 2 else [])
    return complex(np.linalg.lstsq(np.column_stack(cols), np.asarray(z), rcond=None)[0][1])


def timed_psi(beta, V, zt):
    """(ML amplitude z, psi) of one arm; mean mu + z e^{iwt}. With two time levels the means are saturated and psi is the
    within-time-point scatter beyond the Poisson covariance (same ML routine, no change)."""
    t = np.asarray(zt, float) * W
    _, _, z, psi, _ = factor._known_f(beta, V, np.column_stack([np.cos(t), np.sin(t)]), np.zeros(len(t), int))
    return complex(z[0]), float(psi[0])


def arm_stats(beta, V, zt):
    return abs(amp(beta[:, 0] + 1j * beta[:, 1], zt)) ** 2, timed_psi(beta, V, zt)[1]


def kappa(dpsi, dA, comp=0.0):
    num, den = np.nanmedian(np.asarray(dpsi, float) - comp), np.nanmedian(dA)
    return float(num / den) if np.isfinite(num) and np.isfinite(den) and den != 0 else np.nan


def kappas(dp, dA, grp, comp=0.0):
    """(kappa_O, kappa_I) from per-unit dpsi, dA (arrays over units), group labels `grp`, companion shift `comp`."""
    dp, dA, comp = np.asarray(dp, float), np.asarray(dA, float), np.broadcast_to(comp, dp.shape)
    o, i = grp == "O", grp == "I"
    return kappa(dp[o], dA[o], comp[o]), kappa(dp[i], dA[i], comp[i])


# ---- blocks (tissue of a dataset) and units ------------------------------------------------------------------------

def load_blocks():
    """{(dataset, tissue): block} with states of all libraries of the tissue, and {dataset: gene QC table}."""
    blocks, qcs = {}, {}
    for name, (tissues, acol, tcol, qcc) in DATASETS.items():
        pb, donors = data.load(name)
        if name == "vipko_gse296764":
            donors = donors.assign(arm=donors.genotype + "_" + donors.feeding)
        qc = QC.gene_qc(pb if qcc is None else pb[pb.compartment.isin(qcc)], name)
        qcs[name] = qc
        wf = QC.waveform_without(sorted(qc[~qc.kept].gene.unique()))
        for t in tissues:
            x = pb[pb.compartment == t].set_index("donor")
            beta, V = states_cov(x, wf, x.depth.to_numpy(float))
            d = donors.reindex(x.index)
            arm, zt = d[acol].astype(str).to_numpy(), d[tcol].to_numpy(float)
            if name != "vipko_gse296764":  # states equal A297's (same projector, same genes)
                ref, _ = A297.liver_states(name, t, "arm")
                assert np.allclose(ref.z.to_numpy(), beta[:, 0] + 1j * beta[:, 1]) and (ref.index == x.index).all()
            blocks[(name, t)] = dict(dataset=name, tissue=t, beta=beta, V=V, zt=zt, arm=arm, wf=wf, depth=x.depth.to_numpy(float),
                                     names=x.index.astype(str).to_numpy())
    return blocks, qcs


def make_units(blocks):
    units = []
    for name, t, grp, ctl, con, feed in UNITS:
        b = blocks[(name, t)]
        u = dict(dataset=name, tissue=t, group=grp, ctl=ctl, con=con, feeding=feed, block=(name, t))
        for w in ("ctl", "con"):
            r = np.flatnonzero(b["arm"] == u[w])
            u[w + "_rows"], u[w + "_names"], u[w + "_zt"], u[w + "_V"] = r, b["names"][r], b["zt"][r], b["V"][r]
            u[w + "_beta"] = b["beta"][r]
        units.append(u)
    return units


def precondition(u, rng):
    """A395 on the control arm: timed |a| and the q95 of N_PERM time-label permutations (reads the control arm only)."""
    z = u["ctl_beta"][:, 0] + 1j * u["ctl_beta"][:, 1]
    obs = abs(amp(z, u["ctl_zt"]))
    null = np.array([abs(amp(z, rng.permutation(u["ctl_zt"]))) for _ in range(N_PERM)])
    q95 = float(np.quantile(null, 0.95))
    return obs, q95, bool(obs > q95)


# ---- bootstrap (joint over units) -----------------------------------------------------------------------------------

def groups_of(units):
    """Per dataset: {(arm, time): donor names present in every tissue block of the dataset}, plus name -> position maps."""
    gs = {}
    for ds in DATASETS:
        us = [u for u in units if u["dataset"] == ds]
        if not us:
            continue
        g = {}
        for u in us:
            for w in ("ctl", "con"):
                for z in np.unique(u[w + "_zt"]):
                    nm = set(u[w + "_names"][u[w + "_zt"] == z])
                    g[(u[w], z)] = nm if (u[w], z) not in g else g[(u[w], z)] & nm
        gs[ds] = {k: np.array(sorted(v)) for k, v in g.items()}
    return gs


def boot_draw(units, gs, rng):
    """One bootstrap draw: libraries resampled within (arm, time point) by donor name, the same donors in every tissue of a
    dataset (Manella: the same mouse in the three tissues; VIP-KO: liver and heart of the same mouse). Returns per unit
    (dpsi, dA) from the units' current betas."""
    samp = {ds: {k: rng.choice(v, len(v)) for k, v in g.items()} for ds, g in gs.items()}
    dp, dA = np.full(len(units), np.nan), np.full(len(units), np.nan)
    for j, u in enumerate(units):
        st = []
        for w in ("ctl", "con"):
            pos = {n: i for i, n in enumerate(u[w + "_names"])}
            idx = np.array([pos[n] for z in np.unique(u[w + "_zt"]) for n in samp[u["dataset"]][(u[w], z)] if n in pos], int)
            st.append(arm_stats(u[w + "_beta"][idx], u[w + "_V"][idx], u[w + "_zt"][idx]))
        dA[j], dp[j] = lratio(st[0][0], st[1][0]), lratio(st[0][1], st[1][1])
    return dp, dA


# ---- genome-wide tables for the companions ----------------------------------------------------------------------------

def _tidy_depth(x: pd.DataFrame, depth: pd.Series):
    x = x.copy()
    x.index = [canonical(g) for g in x.index]
    return x.groupby(level=0).sum().T, depth.reindex(x.columns)


def genome_greenwell():
    from common.datasets.greenwell2019 import build as G
    raw = G.read_raw()
    depth = raw.sum()
    sym = mouse_ids.to_symbols(raw.index)
    ok = sym.notna().to_numpy()
    x = raw[ok]
    x.index = sym[ok].to_numpy()
    donor = G.meta_of(raw.columns).donor
    x.columns, depth.index = donor.reindex(x.columns).to_numpy(), donor.reindex(depth.index).to_numpy()
    return {"liver": _tidy_depth(x, depth)}


def genome_vipko():
    from common.datasets.vipko_gse296764 import build as V
    import re
    x = pd.read_excel(V.SRC / "GSE296764_SCN_Liver_Heart_count.xlsx", engine="openpyxl")
    x = x.loc[:, ~x.columns.duplicated()].set_index("name")
    cols = [c for c in x.columns if re.fullmatch(r"S\d+(gan|xin|SCN)", c)]
    cnt = x[cols]
    depth = cnt.sum()
    sym = mouse_ids.to_symbols(cnt.index)
    ok = sym.notna().to_numpy()
    out = {}
    for t, suf in (("liver", "gan"), ("heart", "xin")):
        c = [k for k in cols if k.endswith(suf)]
        xt = cnt.loc[ok, c]
        xt.index = sym[ok].to_numpy()
        names = [k[:-len(suf)] for k in c]
        d = depth[c].copy()
        xt.columns, d.index = names, names
        out[t] = _tidy_depth(xt, d)
    return out


def companions(blocks, units, t0, n_sets):
    """dpsi and dA of `n_sets` expression-matched sets per dataset for every unit (rows: dataset, tissue, units' arms)."""
    gen = {"manella2021": A415.genome_manella, "greenwell2019": genome_greenwell, "vipko_gse296764": genome_vipko}
    rows = []
    for di, name in enumerate(DATASETS):
        us = [u for u in units if u["dataset"] == name]
        if not us:
            continue
        g = gen[name]()
        for (ds, t), b in blocks.items():  # the genome tables agree with the pseudobulk depths and libraries
            if ds == name and t in g:
                x, d = g[t]
                x, d = x.reindex(b["names"]), d.reindex(b["names"])
                assert not x.isna().all(axis=1).any() and np.allclose(d.to_numpy(float), b["depth"]), (name, t)
                g[t] = (x, d)
        allx = pd.concat([x for x, _ in g.values()])
        alld = pd.concat([d for _, d in g.values()]).reindex(allx.index)
        lm = pd.Series(np.log((allx.to_numpy(np.float32) / alld.to_numpy(float)[:, None]).mean(0) + 1e-12), index=allx.columns)
        lm = lm[allx.mean(0).to_numpy() >= 1]
        miss = [c for c in WF.genes if c not in lm.index]
        if miss:
            lm = pd.concat([lm, pd.Series(float(lm.median()), index=miss)])
        sets = pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, di]), WF, "matched")
        print(f"companions {name}: {len(lm)} matching genes, {len(sets)} sets ({time.time() - t0:.0f}s)", flush=True)
        fails = {}
        for tk, (xt, _) in g.items():
            mt = xt.median()
            for w in sets:
                fails.setdefault(w.name, set()).update(c for c in w.genes if mt.get(c, 0) < QC.EXPR_MIN)
        for w in sets:
            wk = w.drop(sorted(fails[w.name])) if fails[w.name] else w
            cache = {}
            for u in us:
                try:
                    if u["tissue"] not in cache:
                        x, d = g[u["tissue"]]
                        cache[u["tissue"]] = states_cov(x[list(wk.genes)], wk, d.to_numpy(float))
                    beta, V = cache[u["tissue"]]
                    b = blocks[u["block"]]
                    st = [arm_stats(beta[u[k + "_rows"]], V[u[k + "_rows"]], b["zt"][u[k + "_rows"]]) for k in ("ctl", "con")]
                    rows.append(dict(dataset=name, tissue=u["tissue"], ctl=u["ctl"], con=u["con"], set=w.name, n_genes=len(wk.genes),
                                     dpsi=lratio(st[0][1], st[1][1]), dA=lratio(st[0][0], st[1][0]), error=""))
                except Exception as e:
                    rows.append(dict(dataset=name, tissue=u["tissue"], ctl=u["ctl"], con=u["con"], set=w.name, n_genes=0,
                                     dpsi=np.nan, dA=np.nan, error=f"{type(e).__name__}: {e}"[:120]))
        print(f"companions {name} done ({time.time() - t0:.0f}s)", flush=True)
    c = pd.DataFrame(rows)
    c["unit"] = c.dataset + ":" + c.tissue + ":" + c.ctl + ">" + c.con
    s = c.groupby("unit").agg(sets=("set", "nunique"), finite=("dpsi", lambda v: int(np.isfinite(v).sum())),
                              failed=("error", lambda v: int((v != "").sum())),
                              min_genes=("n_genes", lambda v: int(v[v > 0].min()) if (v > 0).any() else 0)).reset_index()
    return c, s


def comp_matrix(c, units):
    m = c.pivot(index="set", columns="unit", values="dpsi")
    return m[[ulabel(u) for u in units]].to_numpy(float)


# ---- G417.1 simulation ------------------------------------------------------------------------------------------------

def simulate(units, inp, cm, grp, truth, rep, gs):
    """One replicate: simulate both arms of every unit at the observed design (amplitudes, control psi from the data);
    truth S: O perturbed psi scales with |a|^2, I constant; truth N: both constant. Returns (kappa_O, kappa_I, lo90, hi90 of
    kappa_O - kappa_I over N_BOOT_SIM bootstraps, median dA)."""
    rng = np.random.default_rng([SEED, ord(truth), rep])
    full = np.nanmedian(cm, axis=0)
    sim = []
    for u, (z0, p0, z1) in zip(units, inp):
        psi1 = p0 * abs(z1) ** 2 / abs(z0) ** 2 if (truth == "S" and u["group"] == "O") else p0
        v = dict(u)
        for w, z, p in (("ctl", z0, p0), ("con", z1, psi1)):
            r = u[w + "_rows"]
            v[w + "_beta"] = sim_arm(dict(zt=u[w + "_zt"], V=u[w + "_V"]), np.arange(len(r)), z, p, rng)
        sim.append(v)
    dp, dA = np.full(len(sim), np.nan), np.full(len(sim), np.nan)
    for j, v in enumerate(sim):
        s0, s1 = (arm_stats(v[w + "_beta"], v[w + "_V"], v[w + "_zt"]) for w in ("ctl", "con"))
        dA[j], dp[j] = lratio(s0[0], s1[0]), lratio(s0[1], s1[1])
    ko, ki = kappas(dp, dA, grp)
    diffs = []
    for _ in range(N_BOOT_SIM):
        bp, ba = boot_draw(sim, gs, rng)
        sub = rng.choice(len(cm), min(N_SUB, len(cm)), replace=False)
        a, b = kappas(bp, ba, grp, np.nanmedian(cm[sub], axis=0) - full)
        diffs.append(a - b)
    lo, hi = np.nanquantile(diffs, [0.05, 0.95])
    return dict(truth=truth, rep=rep, kappa_O=ko, kappa_I=ki, diff=ko - ki, lo90=float(lo), hi90=float(hi),
                above0=bool(lo > 0), median_dA_O=float(np.nanmedian(dA[grp == "O"])), median_dA_I=float(np.nanmedian(dA[grp == "I"])))


_CTX = {}


def _job(args):
    return simulate(*_CTX["a"], *args)


def power(a):
    limit_threads(WORKERS)
    t0 = time.time()
    sfx = ("_allunits" if a.keep_all else "") + ("_check" if a.check else "")
    OUT.mkdir(parents=True, exist_ok=True)
    blocks, qcs = load_blocks()
    allunits = make_units(blocks)
    rng = np.random.default_rng([SEED, 0])
    rows = []
    for u in allunits:
        obs, q95, ok = precondition(u, rng)
        rows.append(dict(unit=ulabel(u), group=u["group"], feeding=u["feeding"], n_ctl=len(u["ctl_rows"]), n_con=len(u["con_rows"]),
                         n_ctl_times=len(np.unique(u["ctl_zt"])), ctl_abs_a=obs, ctl_q95=q95, a395_pass=ok))
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / f"power_units{sfx}.csv", index=False)
    pd.set_option("display.width", 220)
    print(tab.round(4).to_string(index=False), flush=True)
    keep = np.ones(len(allunits), bool) if a.keep_all else tab.a395_pass.to_numpy()
    units = [u for u, k in zip(allunits, keep) if k]
    grp = np.array([u["group"] for u in units])
    print(f"units kept {len(units)}/{len(allunits)} (O {int((grp == 'O').sum())}, I {int((grp == 'I').sum())}); dropped: "
          + ", ".join(ulabel(u) for u, k in zip(allunits, keep) if not k), flush=True)
    if not (grp == "O").any() or not (grp == "I").any():
        pd.Series(dict(units=len(units), G417_1_decidable=False, reason="an arm has no unit")).to_csv(OUT / f"power_summary{sfx}.csv", header=["value"])
        print("an arm has no unit: not decidable", flush=True)
        return
    n_sets = 8 if a.check else N_SETS
    n_rep = 4 if a.check else N_REP
    c, s = companions(blocks, units, t0, n_sets)
    c.to_csv(OUT / f"companions{sfx}.csv", index=False)
    s.to_csv(OUT / f"companions_summary{sfx}.csv", index=False)
    print(s.to_string(index=False), flush=True)
    cm = comp_matrix(c, units)
    inp = []
    for u in units:
        z0, p0 = timed_psi(u["ctl_beta"], u["ctl_V"], u["ctl_zt"])
        z1, _ = timed_psi(u["con_beta"], u["con_V"], u["con_zt"])
        inp.append((z0, p0, z1))
    gs = groups_of(units)
    _CTX["a"] = (units, inp, cm, grp)
    jobs = [(t, r, gs) for t in "SN" for r in range(n_rep)]
    with mp.get_context("fork").Pool(WORKERS) as pool:
        res = []
        for i, r in enumerate(pool.imap(_job, jobs)):
            res.append(r)
            if (i + 1) % 10 == 0:
                print(f"replicates {i + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    sim = pd.DataFrame(res)
    sim.to_csv(OUT / f"power_sim{sfx}.csv", index=False)
    S, N = sim[sim.truth == "S"], sim[sim.truth == "N"]
    summ = dict(units=len(units), units_O=int((grp == "O").sum()), units_I=int((grp == "I").sum()), n_sets=n_sets, n_rep=n_rep,
                n_boot_sim=N_BOOT_SIM, S_share_above0=float(S.above0.mean()), N_share_above0=float(N.above0.mean()),
                G417_1_S_pass=bool(S.above0.mean() >= 0.8), G417_1_N_pass=bool(N.above0.mean() <= 0.10),
                S_kappa_O_median=float(S.kappa_O.median()), S_kappa_O_sd=float(S.kappa_O.std(ddof=1)),
                S_kappa_I_median=float(S.kappa_I.median()), S_kappa_I_sd=float(S.kappa_I.std(ddof=1)),
                N_kappa_O_median=float(N.kappa_O.median()), N_kappa_O_sd=float(N.kappa_O.std(ddof=1)),
                N_kappa_I_median=float(N.kappa_I.median()), N_kappa_I_sd=float(N.kappa_I.std(ddof=1)),
                S_diff_median=float(S["diff"].median()), n_nan_kappa=int(sim[["kappa_O", "kappa_I"]].isna().sum().sum()))
    summ["G417_1_decidable"] = bool(summ["G417_1_S_pass"] and summ["G417_1_N_pass"])
    pd.Series(summ).to_csv(OUT / f"power_summary{sfx}.csv", header=["value"])
    print(pd.Series(summ).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


# ---- score --------------------------------------------------------------------------------------------------------------

def holm(p):
    order = np.argsort(p)
    adj, run = np.empty(len(p)), 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - r) * p[i]))
        adj[i] = run
    return adj


def score(a):
    limit_threads(WORKERS)
    t0 = time.time()
    sfx = "_check" if a.check else ""
    ps = pd.read_csv(OUT / "power_summary.csv", index_col=0).value
    decidable = str(ps["G417_1_decidable"]) == "True"
    pu = pd.read_csv(OUT / "power_units.csv")
    blocks, _ = load_blocks()
    units = [u for u in make_units(blocks) if bool(pu.set_index("unit").a395_pass[ulabel(u)])]
    c = pd.read_csv(OUT / "companions.csv")
    cm = comp_matrix(c, units)
    cmed = np.nanmedian(cm, axis=0)
    grp = np.array([u["group"] for u in units])
    vip = np.array([u["dataset"] == "vipko_gse296764" for u in units])
    feed = np.array([u["feeding"] for u in units])
    dp, dA = np.full(len(units), np.nan), np.full(len(units), np.nan)
    for j, u in enumerate(units):
        s0, s1 = (arm_stats(u[w + "_beta"], u[w + "_V"], u[w + "_zt"]) for w in ("ctl", "con"))
        dA[j], dp[j] = lratio(s0[0], s1[0]), lratio(s0[1], s1[1])
    gs = groups_of(units)
    rng = np.random.default_rng([SEED, 99])
    n_boot = 200 if a.check else N_BOOT
    BD, BA, KS = [], [], []
    for b in range(n_boot):
        bp, ba = boot_draw(units, gs, rng)
        BD.append(bp), BA.append(ba)
        KS.append(rng.choice(len(cm), min(N_SUB, len(cm)), replace=False))
        if (b + 1) % 250 == 0:
            print(f"bootstrap {b + 1}/{n_boot} ({time.time() - t0:.0f}s)", flush=True)
    BD, BA = np.array(BD), np.array(BA)
    cs_ = np.array([np.nanmedian(cm[s], axis=0) for s in KS])

    def boot(mask_o, mask_i):
        ko, ki = [], []
        for j in range(n_boot):
            ko.append(kappa(BD[j][mask_o], BA[j][mask_o], cs_[j][mask_o]))
            ki.append(kappa(BD[j][mask_i], BA[j][mask_i], cs_[j][mask_i]))
        return np.array(ko), np.array(ki)

    rows = []
    for lab, mo, mi in (("all", grp == "O", grp == "I"), ("I_without_VIPKO", grp == "O", (grp == "I") & ~vip),
                        *((f"feeding_{f}", (grp == "O") & (feed == f), (grp == "I") & (feed == f)) for f in ("AL", "DF", "TRF"))):
        if not mo.any() or not mi.any():
            continue
        ko_obs, ki_obs = kappa(dp[mo], dA[mo], cmed[mo]), kappa(dp[mi], dA[mi], cmed[mi])
        ko, ki = boot(mo, mi)
        d = ko - ki
        q = lambda v: np.nanquantile(v, [0.05, 0.95])  # noqa: E731
        rows.append(dict(subset=lab, units_O=int(mo.sum()), units_I=int(mi.sum()), kappa_O=ko_obs, kappa_O_lo=q(ko)[0], kappa_O_hi=q(ko)[1],
                         kappa_I=ki_obs, kappa_I_lo=q(ki)[0], kappa_I_hi=q(ki)[1], diff=ko_obs - ki_obs, diff_lo=q(d)[0], diff_hi=q(d)[1],
                         p_diff_le0=float(np.nanmean(d <= 0)), p_O_le0=float(np.nanmean(ko <= 0)),
                         p_I_outside=float(max(np.nanmean(ki >= 0.5), np.nanmean(ki <= -0.5))), n_boot_ok=int(np.isfinite(d).sum())))
    out = pd.DataFrame(rows)
    r = out[out.subset == "all"].iloc[0]
    h = holm(np.array([r.p_O_le0, r.p_I_outside]))
    out["H417_1"] = ["pass" if (r.diff_lo > 0) else "fail" if s == "all" else "" for s in out.subset]
    out["H417_2_I_inside"] = [bool(h[1] < 0.05) if s == "all" else None for s in out.subset]
    out["H417_3_O_above0"] = [bool(h[0] < 0.05) if s == "all" else None for s in out.subset]
    out["status"] = ["registered" if decidable and s == "all" else "descriptive (G417.1 not met)" if s == "all" else "descriptive" for s in out.subset]
    tab = pd.DataFrame([dict(unit=ulabel(u), group=u["group"], n_ctl=len(u["ctl_rows"]), n_con=len(u["con_rows"]), dA=dA[j], dpsi=dp[j],
                             comp_dpsi=cmed[j], dpsi_star=dp[j] - cmed[j]) for j, u in enumerate(units)])
    tab.to_csv(OUT / f"score_units{sfx}.csv", index=False)
    out.to_csv(OUT / f"score_kappa{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(tab.round(3).to_string(index=False), out.round(3).T.to_string(), sep="\n")


# ---- descriptive units: loadability only --------------------------------------------------------------------------------

def loadability():
    rows = []
    D = DATA
    try:
        s = pd.read_csv(D / "vollmers2009_GSE13093/samples.csv")
        e = pd.read_csv(D / "vollmers2009_GSE13093/expression_clock.csv", nrows=5)
        rows.append(dict(unit="Vollmers 2009 Cry DKO liver", loadable=True, detail=f"samples.csv {s.shape} cols {list(s.columns)[:8]}; "
                         f"expression_clock.csv cols {e.shape[1]} (microarray, no counts: companions only)"))
    except Exception as ex:
        rows.append(dict(unit="Vollmers 2009 Cry DKO liver", loadable=False, detail=f"{type(ex).__name__}: {ex}"[:200]))
    try:
        import anndata
        ad = anndata.read_h5ad(D / "aorta_jetlag_auerbach2026/adata_qc_filtered_compressed.h5ad", backed="r")
        rows.append(dict(unit="Auerbach 2026 aorta Bmal1 KO", loadable=True, detail=f"h5ad {ad.shape}; obs columns {list(ad.obs.columns)[:12]}"))
    except Exception as ex:
        rows.append(dict(unit="Auerbach 2026 aorta Bmal1 KO", loadable=False, detail=f"{type(ex).__name__}: {ex}"[:200]))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "descriptive_loadability.csv", index=False)
    print(t.to_string(index=False), flush=True)


def main():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--keep-all", action="store_true", help="power sensitivity: simulate all units, ignoring the A395 precondition")
    a = ap.parse_args()
    if a.stage == "power":
        loadability()
    (power if a.stage == "power" else score)(a)


if __name__ == "__main__":
    main()
