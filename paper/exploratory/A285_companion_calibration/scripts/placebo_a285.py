"""A285: placebo calibration and calibrated restatement of the companion tests of the bulk and tumour results.

Registered in paper/docs/registered_experiments.md (A285). Per test: the clock's statistic and its 200 + 200 companions
exactly as scored (states computed once by the scored script's functions); 2000 placebo draws re-pair every set
identically (pairing unit permuted within the scored strata; seeds [285, test, d]). Draw 0 (no permutation) must
reproduce the scored statistic and companion q95.

  V285  placebo rejection rate of the scored rule (clock > companion q95) <= 0.0646, per test x family
  H285  p_cal = (1 + #{placebo pi >= pi}) / (1 + draws), pi = the clock's percentile among its companions

    python placebo_a285.py --test a251|a255|a138|a227 [--draws 2000] [--jobs 4]
    python placebo_a285.py --test ... --families matched,responsive,responsive_broad --out a320   # A320
Outputs (../output/[<out>/]<test>/): repro.csv, draws.npz, tests.csv (A320 also overlap.csv). Defaults reproduce A285.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync.pairwise import calibrated_p, percentile  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, BAND = 285, 0.0646
TEST_CODE = {"a251": 1, "a255": 2, "a138": 3, "a227": 4}
FAMILIES = ["matched", "responsive"]  # set by --families; seeds per family below (scored ones unchanged)
TAG = {"a251": {"matched": 3, "responsive": 4, "responsive_broad": 5},
       "a255": {"matched": 3, "responsive": 4, "responsive_broad": 5},
       "a138": {"matched": 44, "responsive": 45, "responsive_broad": 46},
       "a227": {"matched": [44], "responsive": [44], "responsive_broad": [45]}}
OVERLAP = []  # mean pairwise Jaccard overlap of the gene sets per test x family (where computed)


def jaccard(wfs):
    g = [set(w.genes) for w in wfs]
    j = [len(a & b) / len(a | b) for i, a in enumerate(g) for b in g[:i]]
    return float(np.mean(j)) if j else np.nan


def slices(n_sets, start=1):
    return {f: slice(start + i * n_sets, start + (i + 1) * n_sets) for i, f in enumerate(FAMILIES)}


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def log(*a):
    print(*a, flush=True)


# ---- generic: per-unit vectors (sets x units), NaN where a unit lacks the compartment -----------------------------

def pooled_r(Z: dict, pairs: list, perm: dict | None = None) -> np.ndarray:
    """Weighted mean over pairs of |r| (weights = units with both), for every set at once (A251 / A255 `pooled`)."""
    if perm is not None:
        Z = {c: v[:, perm[c]] for c, v in Z.items()}
    num, den = 0.0, 0.0
    for a, b in pairs:
        ok = ~(np.isnan(Z[a][0]) | np.isnan(Z[b][0]))
        A, B = Z[a][:, ok], Z[b][:, ok]
        r = np.abs(np.sum(A * np.conj(B), 1)) / np.sqrt(np.sum(np.abs(A) ** 2, 1) * np.sum(np.abs(B) ** 2, 1))
        num, den = num + ok.sum() * r, den + ok.sum()
    return num / den


def readout(name, obs, null, families, rule_families):
    """Rows per family: scored rule under placebo (V285) and the calibrated restatement (H285)."""
    rows = []
    for kind, sl in families.items():
        q95 = np.nanquantile(null[:, sl], 0.95, axis=1)
        rate = float(np.mean(null[:, 0] > q95))
        o, nl = percentile(obs, sl)[0], percentile(null, sl)
        rows.append(dict(test=name, family=kind, clock=float(obs[0]), companion_median=float(np.nanmedian(obs[sl])),
                         companion_q95=float(np.nanquantile(obs[sl], 0.95)), scored_rule_pass=bool(obs[0] > np.nanquantile(obs[sl], 0.95)),
                         placebo_rate=rate, V285_pass=rate <= BAND, pi=float(o), p_cal=float(calibrated_p(o, nl)[0]),
                         placebo_rate_pi_1=float(np.mean(nl == 1.0)), in_rule=kind in rule_families))
    return rows


# ---- A251 -----------------------------------------------------------------------------------------------------

def a251():
    M = _load("a285_a251", PAPER / "exploratory/A251_skin_sites/scripts/compute_a251.py")
    from clocksync import paired as P
    data = M.load()
    ids = sorted(set(data["epidermis"][1].query("arm != 'timed_men'").subject_id)
                 | set(data["dermis"][1].query("arm != 'timed_men'").subject_id))
    lm = (sum(data[l][0].loc[data[l][1].arm != "timed_men"].mean() for l in M.LAYERS) / 2 * np.log(2)).dropna()
    lm = lm[~lm.index.duplicated()]
    wfs, fam = [M.WF], {}
    for kind in FAMILIES:
        n0 = len(wfs)
        for w in P.pseudo_waveforms(lm, M.N_SETS, np.random.default_rng([M.SEED, TAG["a251"][kind]]), M.WF, kind):
            if M.women_states(data, w, ids) is not None:
                wfs.append(w)
        fam[kind] = slice(n0, len(wfs))
        OVERLAP.append(dict(test="a251", family=kind, jaccard=jaccard(wfs[n0:])))
    states = [M.women_states(data, w, ids) for w in wfs]
    Z = {c: np.array([s[c] for s in states]) for c in M.COMPS}
    tests = {k.split()[0]: p for k, p in M.families().items()}
    n = len(ids)

    def draw(d):
        if d is None:
            return {k: pooled_r(Z, p) for k, p in tests.items()}
        rng = np.random.default_rng([SEED, TEST_CODE["a251"], d])
        perm = {c: rng.permutation(n) for c in Z}
        return {k: pooled_r(Z, p, perm) for k, p in tests.items()}

    return tests, fam, draw, {k: tuple(FAMILIES) for k in tests}


# ---- A255 -----------------------------------------------------------------------------------------------------

def a255():
    M = _load("a285_a255", PAPER / "exploratory/A255_living_tissues/scripts/compute_a255.py")
    gates = pd.read_csv(PAPER / "exploratory/A255_living_tissues/output/gates.csv")
    mu, _ = M.load_muther(False)
    surg = {a: M.load_surg(a, False) for a in M.SURG if not (gates[gates.cohort == a].p >= 0.05).any()}

    def stack(c, tissues):
        wfs = [c.wf]
        for f in FAMILIES:
            ps = c.pseudo(f, M.N_SETS, TAG["a255"][f])
            OVERLAP.append(dict(test=f"a255 {c.name}", family=f, jaccard=jaccard(ps)))
            wfs += ps
        zs = [c.Z(w) for w in wfs]
        return {k: np.array([z[k] for z in zs]) for k in tissues}

    Zm = stack(mu, ("fat", "skin"))
    sp = {a: tuple(c.t)[:2] for a, c in surg.items()}
    Zs = {a: stack(c, sp[a]) for a, c in surg.items()}
    w = np.array([M.r_pair(Zs[a][sp[a][0]][0], Zs[a][sp[a][1]][0])[2] for a in surg], float)
    w = w / w.sum()  # scored weights (persons with both tissues on the real pairing), fixed under placebo
    fam = slices(M.N_SETS)
    names = list(surg)  # plain names only: the draw closure is shipped to joblib workers
    log("surgical cohorts", names, "weights", w.round(3))
    del surg, mu

    def draw(d):
        rng = None if d is None else np.random.default_rng([SEED, TEST_CODE["a255"], d])
        perm = lambda Z: None if rng is None else {k: rng.permutation(v.shape[1]) for k, v in Z.items()}  # noqa: E731
        out = {"muther fat-skin": pooled_r(Zm, [("fat", "skin")], perm(Zm))}
        out["surgical cohorts pooled"] = sum(wi * pooled_r(Zs[a], [sp[a]], perm(Zs[a])) for wi, a in zip(w, names))
        return out

    return None, fam, draw, {k: tuple(FAMILIES) for k in ("muther fat-skin", "surgical cohorts pooled")}


# ---- A138 -----------------------------------------------------------------------------------------------------

def _a138_cancer(c, i, kind_seeds):
    """Centred tumour / normal states and corrected variances of the clock and every companion set of one cancer."""
    sys.path.insert(0, str(PAPER / "exploratory/A138_cancer_estimability/scripts"))
    import q3_lib as q
    from clocksync.config import limit_threads
    from clocksync.paired import cohort_states, pseudo_waveforms
    limit_threads(1)
    samples = q.sh.load_tcga(None, 0)
    case, _ = q.sh.pair_frames(samples, c)
    ct, cn = q.sh.tcga_full_counts(c, case.index)
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    out = {}
    for kind, tag in kind_seeds:
        out[kind] = []
        for w in pseudo_waveforms(lm, 200, np.random.default_rng([138, tag, 0, i]), q.WF, kind=kind):
            g = list(w.genes)
            bt, vt, bn, vn = cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)
            za, zb = q.centred_complex(bt), q.centred_complex(bn)
            out[kind].append((za, zb, q.tr_corrected(za, vt), q.tr_corrected(zb, vn)))
    return c, out


def a138(jobs=4):
    sys.path.insert(0, str(PAPER / "exploratory/A138_cancer_estimability/scripts"))
    import q3_lib as q
    from joblib import Parallel, delayed
    st = q.tcga_states(None, 0)
    keep = [c for c in st if st[c]["gated"] and q.tr_corrected(st[c]["za"], st[c]["vA"]) > 0
            and q.tr_corrected(st[c]["zb"], st[c]["vB"]) > 0]
    res = dict(Parallel(n_jobs=jobs)(delayed(_a138_cancer)(c, i, tuple((f, TAG["a138"][f]) for f in FAMILIES))
                                     for i, c in enumerate(keep)))
    S = {}
    for c in keep:
        u = st[c]
        sets = [(u["za"], u["zb"], q.tr_corrected(u["za"], u["vA"]), q.tr_corrected(u["zb"], u["vB"]))]
        for f in FAMILIES:
            sets += res[c][f]
        S[c] = dict(za=np.array([x[0] for x in sets]), zb=np.array([x[1] for x in sets]),
                    ea=np.array([x[2] for x in sets]), eb=np.array([x[3] for x in sets]))
    w = np.array([st[c]["n"] for c in keep], float)
    fam = slices(200)
    log("A138 cancers", keep)

    def draw(d):
        rng = None if d is None else np.random.default_rng([SEED, TEST_CODE["a138"], d])
        per = []
        for c in keep:
            x = S[c]
            zb = x["zb"] if rng is None else x["zb"][:, rng.permutation(x["zb"].shape[1])]
            n = x["za"].shape[1]
            with np.errstate(invalid="ignore", divide="ignore"):
                r = np.abs(np.sum(x["za"] * np.conj(zb), 1) / (n - 1)) / np.sqrt(x["ea"] * x["eb"])
            per.append(np.where((x["ea"] > 0) & (x["eb"] > 0), r, np.nan))
        per = np.array(per)
        ww = np.where(np.isfinite(per), w[:, None], 0.0)
        return {"H138.4": np.nansum(per * ww, 0) / ww.sum(0)}

    return None, fam, draw, {"H138.4": tuple(FAMILIES)}


# ---- A227 -----------------------------------------------------------------------------------------------------

A227_DIR = PAPER / "exploratory/A227_tme_compartment_pairing"


def _a227_mod():
    return _load("a285_a227", A227_DIR / "scripts/compute_a227.py")


def _a227_unit(study, kind, families=("matched", "responsive")):
    """Stratum-centred tumour / normal states of the clock and the 200 + 200 companion sets (scored seeds; NaN rows for
    a set without estimable states), and the strata, of one cohort x compartment."""
    M = _a227_mod()
    from clocksync.config import limit_threads
    from clocksync.paired import pseudo_waveforms
    limit_threads(1)
    a, b = M.PAIRS[kind]
    u = M.L.pair_states(M.table(study), a, b)
    G = M.genome_pair(study, (a, b), u["pts"])
    ref = G[G.compartment == b]
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth")]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(M.L.q.WF.genes))]
    za, zb = [u["za"]], [u["zb"]]
    for fam in families:
        for w in pseudo_waveforms(lm, 200, np.random.default_rng(M.seed_of(kind, study) + TAG["a227"][fam]), M.L.q.WF, fam):
            p = M.L.pair_states(G, a, b, wf=w, pts=u["pts"])
            za.append(p["za"] if p else np.full(u["n"], np.nan)), zb.append(p["zb"] if p else np.full(u["n"], np.nan))
    return study, kind, dict(za=np.array(za), zb=np.array(zb), stratum=u["stratum"], n=u["n"])


def _raw(za, zb):
    za, zb = za - za.mean(1, keepdims=True), zb - zb.mean(1, keepdims=True)
    return np.abs(np.sum(za * np.conj(zb), 1)) / np.sqrt(np.sum(np.abs(za) ** 2, 1) * np.sum(np.abs(zb) ** 2, 1))


def a227(jobs=2):
    from joblib import Parallel, delayed
    M = _a227_mod()
    kinds = ("myeloid", "T_NK")
    cohorts = M.cohort_list()
    keys = [(s, k) for k in kinds for s in cohorts if M.units_of(s, k)]
    log("A227 units", len(keys))
    res = Parallel(n_jobs=jobs)(delayed(_a227_unit)(s, k, tuple(FAMILIES)) for s, k in keys)
    U = {(s, k): u for s, k, u in res}
    del M
    fam = slices(200)
    names = {"myeloid": "H227.1", "T_NK": "H227.2"}

    def draw(d):
        rng = None if d is None else np.random.default_rng([SEED, TEST_CODE["a227"], d])
        out = {}
        for kind in kinds:
            sel = [x for x in U if x[1] == kind]
            per, w = [], []
            for x in sel:
                u = U[x]
                zb = u["zb"]
                if rng is not None:  # normal-side rows permuted across patients within chemistry stratum
                    idx = np.arange(u["n"])
                    for st in np.unique(u["stratum"]):
                        j = np.flatnonzero(u["stratum"] == st)
                        idx[j] = j[rng.permutation(len(j))]
                    zb = zb[:, idx]
                with np.errstate(invalid="ignore", divide="ignore"):
                    per.append(_raw(u["za"], zb))
                w.append(u["n"])
            per = np.array(per)
            ww = np.where(np.isfinite(per), np.array(w, float)[:, None], 0.0)
            out[names[kind]] = np.nansum(np.nan_to_num(per) * ww, 0) / ww.sum(0)
        return out

    return None, fam, draw, {"H227.1": ("matched",), "H227.2": ("matched",)}


SCORED = {"a251": PAPER / "exploratory/A251_skin_sites/output/score_tests.csv",
          "a255": PAPER / "exploratory/A255_living_tissues/output/score_tests.csv",
          "a138": PAPER / "exploratory/A138_cancer_estimability/output/responsive/h138_4_responsive.csv",
          "a227": PAPER / "exploratory/A227_tme_compartment_pairing/output/tests.csv"}


def scored(test):
    """Scored statistic and companion q95 per test (columns test, r, matched_q95, responsive_q95)."""
    sc = pd.read_csv(SCORED[test])
    if test == "a138":
        get = lambda k: sc[sc.test.str.startswith(k)].iloc[0]  # noqa: E731
        return pd.DataFrame([dict(test="H138.4", r=get("H138.4 registered (perm").statistic,
                                  matched_q95=get("H138.4 registered vs matched").null_q95,
                                  responsive_q95=get("H138.4-R").null_q95)])
    if test == "a227":
        m = sc.set_index("hypothesis")
        r = pd.read_csv(A227_DIR / "output/responsive/tests.csv").set_index("hypothesis")
        return pd.DataFrame([dict(test=h, r=m.statistic[h], matched_q95=m.pseudo_q95[h], responsive_q95=r.pseudo_q95[h])
                             for h in ("H227.1", "H227.2")])
    return sc.rename(columns={sc.columns[0]: "test"})
BUILD = {"a251": a251, "a255": a255, "a138": a138, "a227": a227}


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--test", choices=sorted(TEST_CODE), required=True)
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--families", default="matched,responsive")
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    cli.setup(args)
    FAMILIES[:] = args.families.split(",")
    out = OUT / args.out / args.test if args.out else OUT / args.test
    out.mkdir(parents=True, exist_ok=True)
    tests, fam, draw, rules = BUILD[args.test](args.jobs) if args.test in ("a138", "a227") else BUILD[args.test]()
    obs = draw(None)
    rep = []
    sc = scored(args.test)
    for k, v in obs.items():
        row = sc[sc.test.astype(str).str.startswith(k)].iloc[0]
        rep.append(dict(test=k, clock=float(v[0]), scored=float(row.get("r", np.nan)),
                        **{f"{f}_q95": float(np.nanquantile(v[s], 0.95)) for f, s in fam.items()},
                        **{f"scored_{f}_q95": float(row.get(f"{f}_q95", np.nan)) for f in fam}))
    rep = pd.DataFrame(rep)
    rep.to_csv(out / "repro.csv", index=False)
    log(rep.round(4).to_string())
    assert np.allclose(rep.clock, rep.scored, atol=1e-9), "draw 0 does not reproduce the scored statistic"
    for f in fam:
        if f in ("matched", "responsive"):  # scored families only
            assert np.allclose(rep[f"{f}_q95"], rep[f"scored_{f}_q95"], atol=1e-9), f"draw 0: {f} q95 not reproduced"
    if OVERLAP:
        pd.DataFrame(OVERLAP).to_csv(out / "overlap.csv", index=False)
    from joblib import Parallel, delayed
    res = Parallel(n_jobs=args.jobs)(delayed(draw)(d) for d in range(args.draws))
    null = {k: np.array([r[k] for r in res]) for k in obs}
    np.savez_compressed(out / "draws.npz", **null)
    rows = []
    for k in obs:
        rows += readout(k, obs[k], null[k], fam, rules[k])
    T = pd.DataFrame(rows)
    T.to_csv(out / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    log(T.round(4).to_string())


if __name__ == "__main__":
    main()
