"""A212: can bulk coherence statistics tell a damped / dephased clock (D) from a swamped (S), homogeneously rewired (W) or
heterogeneously rewired (H) one? Plasmode simulation on real tumours (registered in paper/docs/registered_experiments.md, A212).

Bases: B1 = reference arm of A176 H127.1 (SCAN-B, lowest proliferation tertile per PAM50 x purity stratum); B2 = reference arm of
H166.1 (METABRIC ER+, A166 H1 groups). Per stratum, Y (ln 2 x log2 expression of the 14 clock genes) is centred and split into the
clock part C = Y P and the rest E = Y - C, with P the clock-plane projector of the gene-centred waveform design (clocksync
`plane`, the projector the A176 statistics use; C = beta Xc' with beta = Y Xc (Xc'Xc)^-1). Each replicate splits every stratum at
random into halves; the case half receives the mechanism; statistics are the A176 ones (decompose.components: dR, dlog_plane,
dlog_off, ds), case minus reference, mean over strata.
  --stage grid    : 50 replicates per grid value, matched strength per base x mechanism (output/grid.csv, matched.csv)
  --stage clouds  : 200 replicates per mechanism at the matched strength (output/clouds.csv)
  --stage tests   : cross-validated LDA AUC (V212.1-V212.3), observed A176 points against the clouds (output/tests.csv, observed.csv)
Seeds: SeedSequence([crc32("A212"), base, mechanism, stage, grid index, replicate]).
    python compute_a212.py --stage grid|clouds|tests [--jobs 4]"""
import argparse
import os
import sys
import zlib
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
sys.path.insert(0, str(EXP / "A168_endocrine/scripts"))
sys.path.insert(0, str(EXP / "A127_scanb/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold, cross_val_predict  # noqa: E402

import a127_lib as A127  # noqa: E402
import a168_lib as L  # noqa: E402
from clocksync._core.study import plane  # noqa: E402
from clocksync.decompose import components  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["grid", "clouds", "tests"])
ap.add_argument("--jobs", type=int, default=4)
args = ap.parse_args()
OUT = HERE.parent / "output"; OUT.mkdir(exist_ok=True)
SEED = zlib.crc32(b"A212") % (2 ** 31)
HOURS = 24 / (2 * np.pi)
MECHS = ["D", "S", "W", "H"]
BASES = {"B1": "H127.1", "B2": "H166.1"}
N_GRID, N_CLOUD, N_FACTORS, AUC_BAR = 50, 200, 3, 0.9
STATS = ["dR", "dlog_plane", "dlog_off", "ds"]
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


# ---- bases ------------------------------------------------------------------------------------------------------------------------
def base_b1():
    """SCAN-B H127.1 reference arm: per PAM50 x purity stratum (>= 30 tumours), the lowest proliferation tertile (A176 build_all)."""
    m, X = A127.load()
    m["pur_t"] = m.groupby("pam50").purity.transform(lambda s: A127.tertile(s))
    Y = L.Yof(X, A127.GENES); prolif = m.prolif.to_numpy(); out = []
    for _, g in m.groupby(["pam50", "pur_t"]):
        if len(g) >= 30:
            i = np.array([m.index.get_loc(k) for k in g.index])
            out.append(Y[i[A127.tertile(prolif[i]) == 0]])
    return out, A127.WF


def base_b2():
    """METABRIC H166.1 reference arm: ER+ tumours of each A166 H1 group."""
    A166 = L.A166
    mb, Xb = A166.load("metabric"); mb["pos"] = np.arange(len(mb))
    h = A166.build(mb)[0]
    Y = L.Yof(Xb, A166.GENES)
    return [Y[i[h.er[i] == 1]] for i in h.G], L.wf_of(list(A166.GENES))


def split(Y, wf):
    """Centred stratum -> (C, E, beta) with C = Y P (gene-centred plane) and beta the clock state in the gene-centred design."""
    y = Y - Y.mean(0)
    Xc = wf.X - wf.X.mean(0)
    beta = y @ Xc @ np.linalg.inv(Xc.T @ Xc)
    C = y @ plane(wf)
    return C, y - C, beta


def design(wf, shift_h):
    """Gene-centred design with acrophases shifted by shift_h (hours, per gene)."""
    psi = wf.psi + shift_h / HOURS
    X = np.column_stack([wf.amp * np.cos(psi), wf.amp * np.sin(psi)])
    return X - X.mean(0)


def apply(mech, s, C, E, beta, wf, rng, sd_e):
    """Case half after mechanism `mech` at strength s."""
    if mech == "D":
        return s * C + E
    if mech == "S":
        U = rng.normal(0.0, 1.0, (len(C), N_FACTORS)); Ld = rng.normal(0.0, s * sd_e, (C.shape[1], N_FACTORS))
        return C + E + U @ Ld.T
    if mech == "W":
        return beta @ design(wf, rng.normal(0.0, s, C.shape[1])).T + E
    return np.stack([beta[i] @ design(wf, rng.normal(0.0, s, C.shape[1])).T for i in range(len(C))]) + E  # H


def comp(y, wf):
    y = y - y.mean(0)
    c = components(y, np.zeros(y.shape[1]), 1, wf)
    return np.array([c["coherence"], c["log_plane"], c["log_off"]])


def replicate(parts, wf, mech, s, seed, sd_e):
    """One replicate: (dR, dlog_plane, dlog_off, ds), mean over strata."""
    rng = np.random.default_rng(seed); d = []
    if mech == "W":  # one shift vector shared by all case tumours of the replicate
        shift = rng.normal(0.0, s, len(wf.genes)); Xw = design(wf, shift)
    for C, E, beta in parts:
        i = rng.permutation(len(C)); a, b = i[: len(i) // 2], i[len(i) // 2:]
        if mech == "W":
            ya = beta[a] @ Xw.T + E[a]
        else:
            ya = apply(mech, s, C[a], E[a], beta[a], wf, rng, sd_e)
        d.append(comp(ya, wf) - comp(C[b] + E[b], wf))
    d = np.mean(d, 0)
    return np.r_[d, d[1] - d[2]]


def grid_of(mech, sd_e):
    return {"D": np.round(np.arange(1.0, 0.0, -0.1), 2), "S": np.linspace(0.0, 2.0, 10),
            "W": np.arange(0.0, 7.0), "H": np.arange(0.0, 7.0)}[mech]


def run(parts, wf, mech, values, n, tags, sd_e):
    jobs = [(v, r, np.random.SeedSequence([SEED, *tags, k, r])) for k, v in enumerate(values) for r in range(n)]
    res = Parallel(n_jobs=args.jobs)(delayed(replicate)(parts, wf, mech, v, s, sd_e) for v, r, s in jobs)
    return pd.DataFrame([dict(strength=v, rep=r, **dict(zip(STATS, x))) for (v, r, _), x in zip(jobs, res)])


# ---- main -------------------------------------------------------------------------------------------------------------------------
obs = pd.read_csv(EXP / "A176_c3_plane_variance/output/tests.csv").set_index("contrast")
B = {}
for b, f in (("B1", base_b1), ("B2", base_b2)):
    strata, wf = f()
    parts = [split(Y, wf) for Y in strata]
    sd_e = float(np.median(np.concatenate([E for _, E, _ in parts]).std(0)))  # per-gene SD of E (median over genes)
    B[b] = dict(parts=parts, wf=wf, sd_e=sd_e, target=float(obs.loc[BASES[b], "dR"]))
    say(b, "strata", len(parts), "tumours", sum(len(p[0]) for p in parts), "median SD of E", round(sd_e, 4), "target dR", B[b]["target"])

if args.stage == "grid":
    rows, match = [], []
    for bi, (b, r) in enumerate(B.items()):
        for mi, mech in enumerate(MECHS):
            g = grid_of(mech, r["sd_e"])
            t = run(r["parts"], r["wf"], mech, g, N_GRID, [bi, mi, 0], r["sd_e"]).assign(base=b, mechanism=mech)
            rows.append(t)
            mean = t.groupby("strength").dR.mean().reindex(g).to_numpy()
            reach = np.flatnonzero(mean <= r["target"])
            if len(reach):  # first grid value at or below the target, interpolated with the previous one
                k = reach[0]
                s = g[0] if k == 0 else float(np.interp(r["target"], [mean[k], mean[k - 1]], [g[k], g[k - 1]]))
                ok = True
            else:
                s, ok = float(g[np.argmin(mean)]), False
            match.append(dict(base=b, mechanism=mech, target_dR=r["target"], matched_strength=s, reaches_observed=ok,
                              min_mean_dR=float(mean.min()), strength_unit={"D": "rho", "S": "loading SD / SD(E)", "W": "h", "H": "h"}[mech]))
            say(b, mech, "matched", round(s, 3), "reaches", ok)
    pd.concat(rows).to_csv(OUT / "grid.csv", index=False)
    M = pd.DataFrame(match); M.to_csv(OUT / "matched.csv", index=False); print(M.round(4).to_string())

elif args.stage == "clouds":
    M = pd.read_csv(OUT / "matched.csv").set_index(["base", "mechanism"])
    rows = []
    for bi, (b, r) in enumerate(B.items()):
        for mi, mech in enumerate(MECHS):
            s = M.loc[(b, mech), "matched_strength"]
            rows.append(run(r["parts"], r["wf"], mech, [s], N_CLOUD, [bi, mi, 1], r["sd_e"]).assign(base=b, mechanism=mech))
            say(b, mech, "cloud done")
    pd.concat(rows).to_csv(OUT / "clouds.csv", index=False)

else:
    Cl = pd.read_csv(OUT / "clouds.csv"); M = pd.read_csv(OUT / "matched.csv").set_index(["base", "mechanism"])
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)

    def auc(b, m1, m2):
        d = Cl[(Cl.base == b) & Cl.mechanism.isin([m1, m2])]
        x, y = d[["dlog_plane", "dlog_off"]].to_numpy(), (d.mechanism == m1).to_numpy().astype(int)
        p = cross_val_predict(LinearDiscriminantAnalysis(), x, y, cv=cv, method="predict_proba")[:, 1]
        return float(roc_auc_score(y, p))
    rows = []
    for b in B:
        for m1, m2 in (("W", "D"), ("H", "D"), ("S", "D"), ("W", "H"), ("S", "W"), ("S", "H")):
            rows.append(dict(base=b, pair=f"{m1} vs {m2}", auc=auc(b, m1, m2)))
    A = pd.DataFrame(rows).set_index(["base", "pair"]).auc
    tests = []
    for name, b in (("V212.1", "B1"), ("V212.3", "B2")):
        tests.append(dict(test=name, base=b, part="rewiring vs damping", auc_W_D=A[(b, "W vs D")], auc_H_D=A[(b, "H vs D")],
                          passes=bool(min(A[(b, "W vs D")], A[(b, "H vs D")]) >= AUC_BAR)))
    for name, b in (("V212.2", "B1"), ("V212.3", "B2")):
        tests.append(dict(test=name, base=b, part="swamping vs damping", auc_S_D=A[(b, "S vs D")], passes=bool(A[(b, "S vs D")] >= AUC_BAR)))
    T = pd.DataFrame(tests); T.to_csv(OUT / "tests.csv", index=False); A.reset_index().to_csv(OUT / "auc_all.csv", index=False)
    print(T.round(4).to_string()); print(A.round(4).to_string())
    # observed A176 points against the clouds (descriptive; the observed values were seen before registration)
    ob = []
    for b, con in BASES.items():
        o = obs.loc[con, ["dlog_plane", "dlog_off"]].to_numpy(float)
        for mech in MECHS:
            d = Cl[(Cl.base == b) & (Cl.mechanism == mech)][["dlog_plane", "dlog_off"]].to_numpy()
            mu, S = d.mean(0), np.cov(d.T); Si = np.linalg.inv(S)
            md = lambda x: float(np.sqrt((x - mu) @ Si @ (x - mu)))
            ob.append(dict(base=b, contrast=con, mechanism=mech, obs_dlog_plane=o[0], obs_dlog_off=o[1], cloud_dlog_plane=mu[0],
                           cloud_dlog_off=mu[1], cloud_dR=Cl[(Cl.base == b) & (Cl.mechanism == mech)].dR.mean(), mahalanobis=md(o),
                           pct_cloud_closer=100 * float(np.mean([md(x) <= md(o) for x in d])),
                           reaches_observed=bool(M.loc[(b, mech), "reaches_observed"]), matched_strength=M.loc[(b, mech), "matched_strength"]))
    O = pd.DataFrame(ob); O.to_csv(OUT / "observed.csv", index=False); print(O.round(3).to_string())
say("done", args.stage)
