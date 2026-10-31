"""A175 (C1): does the tumour follow its own normal tissue less when it is proliferating or ER-negative?

Registered in paper/docs/registered_experiments.md, section "A175". Reuses the H138.4 / H138.8 library of
`A138_cancer_estimability/scripts` (`q3_lib.tcga_states`, `h138_8.Cancer` / `load_absolute`, `score.pseudo_pooled` rule,
`pseudo_waveforms`) by import. Pairing statistic = H138.4 (|complex correlation| of tumour and own-normal states, states
re-centred within the arm, one projector per cancer); pooled with weights = pairs of the arm.
  H175.1  D_prolif = pooled(low G2M tertile) - pooled(high tertile) > 0; tertiles by rank of the tumour G2M score within
          cancer (A97 definition, z within the paired tumours of the cancer); cancers with >= 10 pairs per arm.
  H175.2  D_ER = pairing(ER+) - pairing(ER-) in BRCA; run only if its MDE <= 0.3 (else "not estimable at this n").
Null: arm labels permuted within cancer (10,000), p = (1 + #{D_null >= D}) / 10,001; Holm over the tests that run.
Power: MDE = q95 + 0.84 SD of that null. The power stage uses labels only and never evaluates the observed D.
Companions (same labels): D on 200 matched pseudo-clock sets (A52 rule) and 200 responsive sets per cancer.

    python run_a175.py --stage power|tests [--jobs 4] [--n-perm 10000] [--n-sets 200]
Outputs (output/): power.csv, power_arms.csv; tests.csv, per_cancer.csv, descriptives.csv, h138_4_reproduction.csv
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parents[1]
A138 = HERE.parent / "A138_cancer_estimability" / "scripts"
sys.path.insert(0, str(A138))
import q3_lib as q  # noqa: E402
from h138_8 import Cancer, load_absolute  # noqa: E402
from common import cli  # noqa: E402
from power import SEED  # noqa: E402  (138)

SCCC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
PAM50 = SCCC / "tcga/pancan_subtypes/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt"
ER_META = HERE.parent / "A166_metabric" / "output" / "tcga_meta.parquet"
MIN_ARM, ER_STOP = 10, 0.3
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:7.1f} min]", *a, flush=True)


# ---- data ----------------------------------------------------------------------------------------------------------

def tertiles(x: np.ndarray) -> np.ndarray:
    """0 low, 1 mid, 2 high by rank within cancer."""
    return pd.qcut(pd.Series(x).rank(method="first"), 3, labels=False).to_numpy()


def g2m_scores(ct: pd.DataFrame) -> np.ndarray:
    """A97 proliferation: mean z of HALLMARK_G2M_CHECKPOINT (clock genes removed) of log2 CPM, z over these tumours."""
    gm = q.sh.hallmark_g2m()
    need = [g for g in gm if g in ct.columns]
    lc = np.log2(ct[need].div(ct.sum(axis=1), axis=0) * 1e6 + 1)
    return q.sh._zmean(lc, gm).to_numpy()


def prepare(st: dict, keep_c: list) -> dict:
    """Per gated cancer: patient ids, G2M score, tertile labels, ABSOLUTE purity, ER and PAM50 labels."""
    samples = q.sh.load_tcga(None, 0)
    absp = load_absolute()
    er = pd.read_parquet(ER_META).drop_duplicates("patient").set_index("patient").er_xml
    pam = pd.read_csv(PAM50, sep="\t", comment="#").set_index("PATIENT_ID").SUBTYPE
    d = {}
    for c in keep_c:
        case, _ = q.sh.pair_frames(samples, c)
        ct, _ = q.sh.tcga_full_counts(c, case.index)
        ct = ct.loc[case.index]
        g = g2m_scores(ct)
        e = dict(cancer=c, patients=case.index.to_numpy(), g2m=g, tert=tertiles(g),
                 purity=case.barcode.str[:15].map(absp).to_numpy(float))
        if c == "BRCA":
            e["er"] = case.index.map(er).map({"Positive": 1, "Negative": 0}).fillna(-1).to_numpy(int)
            e["pam"] = case.index.map(pam).map({"BRCA_LumA": 1, "BRCA_LumB": 1, "BRCA_Basal": 0, "BRCA_Her2": 0}).fillna(-1).to_numpy(int)
        d[c] = e
    return d


def subset(C: Cancer, keep: np.ndarray) -> Cancer:
    S = Cancer.__new__(Cancer)
    S.name, S.n = C.name, int(keep.sum())
    S.a, S.b, S.va, S.vb = C.a[keep], C.b[keep], C.va[keep], C.vb[keep]
    return S


# ---- contrast: pooled(arm a) - pooled(arm b), labels permuted within cancer -----------------------------------------

class Contrast:
    """Labels per cancer (int, -1 = excluded); arms `a` (reference, expected higher) and `b`; cancers with >= MIN_ARM pairs
    in both arms. `null()` uses labels only; `observed()` evaluates the statistic."""

    def __init__(self, cs: dict, labels: dict, a: int, b: int, tag: int, cidx: dict):
        self.a, self.b, self.tag, self.cidx = a, b, tag, cidx
        self.keep = {c: labels[c] >= 0 for c in cs}
        self.lab = {c: labels[c][self.keep[c]] for c in cs}
        self.sub = {c: subset(cs[c], self.keep[c]) for c in cs}
        self.n = {c: (int((self.lab[c] == a).sum()), int((self.lab[c] == b).sum())) for c in cs}
        self.use = [c for c in cs if min(self.n[c]) >= MIN_ARM]

    def _pool(self, s: dict, which: int) -> np.ndarray:
        return q.pool(np.array([s[c][which] for c in self.use]), [self.n[c][which] for c in self.use])

    def null(self, n_perm: int) -> np.ndarray:
        s = {}
        for c in self.use:
            rng = np.random.default_rng([SEED, 175, self.tag, self.cidx[c]])
            perms = np.argsort(rng.random((n_perm, len(self.lab[c]))), axis=1)
            lp = self.lab[c][perms]
            s[c] = (self.sub[c].stat(lp == self.a), self.sub[c].stat(lp == self.b))
        return self._pool(s, 0) - self._pool(s, 1)

    def observed(self, cs: dict | None = None):
        """(per-cancer stats, pooled a, pooled b, D) on the true labels; `cs` replaces the states (companion sets)."""
        s = {}
        for c in self.use:
            C = self.sub[c] if cs is None else subset(cs[c], self.keep[c])
            s[c] = (C.stat((self.lab[c] == self.a)[None, :]), C.stat((self.lab[c] == self.b)[None, :]))
        pa, pb = self._pool(s, 0), self._pool(s, 1)
        return s, float(pa[0]), float(pb[0]), float(pa[0] - pb[0])


def contrasts(cs, d, cidx):
    out = {"H175.1": Contrast(cs, {c: d[c]["tert"] for c in cs}, 0, 2, 1, cidx)}
    br = {"BRCA": cs["BRCA"]}
    out["H175.2"] = Contrast(br, {"BRCA": d["BRCA"]["er"]}, 1, 0, 2, cidx)
    return out


# ---- companions ----------------------------------------------------------------------------------------------------

def companion_cancer(c: str, i: int, n_sets: int, arms: dict) -> dict:
    """D inputs on `n_sets` matched and responsive pseudo-clock sets of cancer `c` (streams of `score.pseudo_pooled` and
    `score_h138_responsive`); `arms`: name -> (labels or None, a, b). Returns stats per set and arm."""
    from clocksync.config import limit_threads
    from clocksync.paired import cohort_states, pseudo_waveforms
    limit_threads(1)
    samples = q.sh.load_tcga(None, 0)
    case, _ = q.sh.pair_frames(samples, c)
    ct, cn = q.sh.tcga_full_counts(c, case.index)
    ct, cn = ct.loc[case.index], cn.loc[case.index]
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    out = {}
    for kind, tag, kw in (("matched", 44, {}), ("responsive", 45, {"kind": "responsive"})):
        ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, tag, 0, i]), q.WF, **kw)
        res = {k: [] for k in arms}
        for w in ws:
            g = list(w.genes)
            bt, vt, bn, vn = cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)
            C = Cancer(c, dict(n=len(case), raw=(bt, vt, bn, vn)))
            for k, (lab, a, b) in arms.items():
                keep = lab >= 0
                S = subset(C, keep)
                res[k].append((float(S.stat((lab[keep] == a)[None, :])[0]), float(S.stat((lab[keep] == b)[None, :])[0])))
        out[kind] = {k: np.array(v) for k, v in res.items()}
    return out


def companion_D(con: Contrast, name: str, per_cancer: list, n_sets: int) -> dict:
    """Pooled D per set index (weights = arm pairs) for each kind."""
    out = {}
    cl = [r for r in per_cancer if r is not None]
    for kind in ("matched", "responsive"):
        pa = q.pool(np.array([r[kind][name][:, 0] for r in cl]), [con.n[c][0] for c in con.use])
        pb = q.pool(np.array([r[kind][name][:, 1] for r in cl]), [con.n[c][1] for c in con.use])
        out[kind] = pa - pb
    return out


# ---- stages --------------------------------------------------------------------------------------------------------

def setup_data():
    st = q.tcga_states(None, 0)
    keep_c = [c for c in st if st[c]["gated"]]
    obs4 = {c: q.corr_stat(st[c]["za"], st[c]["vA"], st[c]["zb"], st[c]["vB"]) for c in keep_c}
    pooled4 = float(q.pool(np.array([[obs4[c]] for c in keep_c]), [st[c]["n"] for c in keep_c])[0])
    n4 = sum(st[c]["n"] for c in keep_c)
    log(f"H138.4 reproduced: pooled {pooled4:.4f}, {len(keep_c)} cancers, {n4} pairs")
    assert abs(pooled4 - 0.3646870339954388) < 5e-4 and n4 == 517 and len(keep_c) == 10, "H138.4 not reproduced"
    cs = {c: Cancer(c, st[c]) for c in keep_c}
    for c in keep_c:
        assert abs(cs[c].stat(np.ones((1, cs[c].n), bool))[0] - obs4[c]) < 1e-9
    cidx = {c: i for i, c in enumerate(keep_c)}
    return st, keep_c, cs, cidx, dict(h138_4_pooled_reproduced=pooled4, n_pairs=n4, n_cancers=len(keep_c))


def power_row(name, con, n_perm):
    null = con.null(n_perm)
    s = q.summarise(null)
    return (dict(hypothesis=name, null_mean=s["null_mean"], null_sd=s["null_sd"], null_q95=s["null_q95"], mde=s["mde"],
                 n_perm=s["n_perm"], cancers=";".join(con.use), n_cancers=len(con.use),
                 n_pairs_arm_a=sum(con.n[c][0] for c in con.use), n_pairs_arm_b=sum(con.n[c][1] for c in con.use),
                 n_pairs_per_arm_per_cancer=";".join(f"{c}:{con.n[c][0]}/{con.n[c][1]}" for c in con.use)),
            null)


def stage_power(args, out):
    st, keep_c, cs, cidx, rep = setup_data()
    d = prepare(st, keep_c)
    cons = contrasts(cs, d, cidx)
    rows, arms = [], []
    for name, con in cons.items():
        r, _ = power_row(name, con, args.n_perm)
        rows.append(r)
        arms += [dict(hypothesis=name, cancer=c, n_arm_a=con.n[c][0], n_arm_b=con.n[c][1], in_pool=c in con.use) for c in con.n]
        log(name, {k: r[k] for k in ("null_mean", "null_sd", "null_q95", "mde", "n_pairs_per_arm_per_cancer")})
    pd.DataFrame(rows).to_csv(out / "power.csv", index=False)
    pd.DataFrame(arms).to_csv(out / "power_arms.csv", index=False)
    log("power stage done")


def pseudo_summary(D, obs):
    x = D[np.isfinite(D)]
    q95 = float(np.quantile(x, 0.95))
    return dict(q95=q95, median=float(np.median(x)), pct=float(100 * np.mean(x <= obs)), exceeds=bool(obs > q95), n=len(x))


def stage_tests(args, out):
    pw = pd.read_csv(out / "power.csv").set_index("hypothesis")
    st, keep_c, cs, cidx, rep = setup_data()
    pd.DataFrame([rep]).to_csv(out / "h138_4_reproduction.csv", index=False)
    d = prepare(st, keep_c)
    cons = contrasts(cs, d, cidx)
    run_er = bool(pw.loc["H175.2", "mde"] <= ER_STOP)
    log(f"H175.2 stop rule: MDE {pw.loc['H175.2', 'mde']:.4f} -> {'run' if run_er else 'not run (not estimable at this n)'}")
    names = ["H175.1"] + (["H175.2"] if run_er else [])
    # companions: sets of every cancer in the pool, labels = observed
    work = {}
    for c in cons["H175.1"].use:
        a = {"H175.1": (d[c]["tert"], 0, 2)}
        if run_er and c == "BRCA":
            a["H175.2"] = (d[c]["er"], 1, 0)
        work[c] = a
    cl = {c: i for i, c in enumerate(keep_c)}
    needed = sorted(set(work) | ({"BRCA"} if run_er else set()), key=cl.get)
    res = Parallel(n_jobs=args.jobs)(delayed(companion_cancer)(c, cl[c], args.n_sets, work.get(c) or {"H175.2": (d[c]["er"], 1, 0)}) for c in needed)
    res = dict(zip(needed, res))
    log("companion sets done")

    rows = []
    for name in names:
        con = cons[name]
        null = con.null(args.n_perm)
        s = q.summarise(null)
        x = null[np.isfinite(null)]
        pw_row = pw.loc[name]
        for k in ("null_mean", "null_sd", "null_q95", "mde"):  # same permutations as the power stage
            assert np.isclose(s[k], pw_row[k], rtol=1e-9, atol=1e-12), f"{name}: null differs from power.csv ({k})"
        _, pa, pb, D = con.observed()
        p = float((1 + np.sum(x >= D)) / (len(x) + 1))
        comp = companion_D(con, name, [res[c] for c in con.use], args.n_sets)
        row = dict(hypothesis=name, D=D, pooled_arm_a=pa, pooled_arm_b=pb, p=p, null_mean=s["null_mean"], null_sd=s["null_sd"],
                   null_q95=s["null_q95"], mde=s["mde"], n_perm=len(x), cancers=";".join(con.use),
                   n_pairs_arm_a=sum(con.n[c][0] for c in con.use), n_pairs_arm_b=sum(con.n[c][1] for c in con.use))
        for kind in ("matched", "responsive"):
            ps = pseudo_summary(comp[kind], D)
            row |= {f"{kind}_q95": ps["q95"], f"{kind}_median": ps["median"], f"{kind}_percentile_of_D": ps["pct"],
                    f"{kind}_D_exceeds_q95": ps["exceeds"], f"{kind}_n_sets": ps["n"]}
        row["clock_specific"] = bool(row["matched_D_exceeds_q95"] and row["responsive_D_exceeds_q95"])
        rows.append(row)
        log(name, {k: row[k] for k in ("D", "p", "mde", "matched_percentile_of_D", "responsive_percentile_of_D")})
    t = pd.DataFrame(rows)
    # Holm over the tests that run
    order = t.p.argsort().to_numpy()
    m, run_max, adj = len(t), 0.0, np.empty(len(t))
    for r, i in enumerate(order):
        run_max = max(run_max, min(1.0, (m - r) * t.p.iloc[i]))
        adj[i] = run_max
    t["p_holm"] = adj
    verd = []
    for _, r in t.iterrows():
        if r.p_holm < 0.05 and r.D > 0:
            v = "supported: lower pairing in the high-proliferation arm" if r.hypothesis == "H175.1" else "supported: lower pairing in ER-"
            v += "; clock-specific" if r.clock_specific else "; not clock-specific (a pseudo-clock companion gives as large a D)"
        elif r.mde > 0.24:
            v = "failed, MDE > 0.24: not resolved at this n"
        else:
            v = "failed, MDE <= 0.24: pairing not lower (coherence lost inside a tumour that still tracks the person)"
        verd.append(v)
    t["verdict"] = verd
    if not run_er:
        t = pd.concat([t, pd.DataFrame([dict(hypothesis="H175.2", mde=pw.loc["H175.2", "mde"],
                                              verdict="not estimable at this n (stop rule: MDE > 0.3), not run")])], ignore_index=True)
    t.to_csv(out / "tests.csv", index=False)

    # descriptives -------------------------------------------------------------------------------------------------
    c1 = cons["H175.1"]
    desc = []
    pcs = []
    for c in keep_c:
        t3 = [float(cs[c].stat((d[c]["tert"] == k)[None, :])[0]) for k in range(3)]
        nn = [int((d[c]["tert"] == k).sum()) for k in range(3)]
        pcs.append(dict(cancer=c, n_pairs=cs[c].n, n_low=nn[0], n_mid=nn[1], n_high=nn[2], stat_low=t3[0], stat_mid=t3[1], stat_high=t3[2],
                        D_low_minus_high=t3[0] - t3[2], in_pool=c in c1.use,
                        spearman_g2m_vs_purity=float(pd.Series(d[c]["g2m"]).corr(pd.Series(d[c]["purity"]), method="spearman"))))
    pd.DataFrame(pcs).to_csv(out / "per_cancer.csv", index=False)
    # mid tertile, pooled over the pool cancers
    mid = q.pool(np.array([[float(cs[c].stat((d[c]["tert"] == 1)[None, :])[0])] for c in c1.use]), [int((d[c]["tert"] == 1).sum()) for c in c1.use])
    desc.append(dict(descriptive="pooled pairing, mid tertile (pool cancers)", value=float(mid[0]), n_pairs=sum(int((d[c]["tert"] == 1).sum()) for c in c1.use)))
    desc.append(dict(descriptive="pooled pairing, low tertile", value=t.loc[0, "pooled_arm_a"], n_pairs=int(t.loc[0, "n_pairs_arm_a"])))
    desc.append(dict(descriptive="pooled pairing, high tertile", value=t.loc[0, "pooled_arm_b"], n_pairs=int(t.loc[0, "n_pairs_arm_b"])))
    # purity-residualised proliferation (ABSOLUTE purity, within cancer)
    lab_r = {}
    for c in keep_c:
        ok = np.isfinite(d[c]["purity"])
        lab = -np.ones(len(ok), int)
        x, y = d[c]["purity"][ok], d[c]["g2m"][ok]
        res_ = y - np.polyval(np.polyfit(x, y, 1), x)
        lab[ok] = tertiles(res_)
        lab_r[c] = lab
    cr = Contrast(cs, lab_r, 0, 2, 3, cidx)
    nr = cr.null(args.n_perm)
    _, pa, pb, Dr = cr.observed()
    x = nr[np.isfinite(nr)]
    desc.append(dict(descriptive="D_prolif, G2M residualised on ABSOLUTE purity", value=Dr, p=float((1 + np.sum(x >= Dr)) / (len(x) + 1)),
                     null_sd=float(x.std(ddof=1)), n_pairs=int(sum((l >= 0).sum() for l in lab_r.values())),
                     note=f"pooled low {pa:.4f}, high {pb:.4f}; cancers {';'.join(cr.use)}; arms {';'.join(f'{c}:{cr.n[c][0]}/{cr.n[c][1]}' for c in cr.use)}"))
    # BRCA PAM50 luminal vs basal + HER2
    cp = Contrast({"BRCA": cs["BRCA"]}, {"BRCA": d["BRCA"]["pam"]}, 1, 0, 4, cidx)
    npm = cp.null(args.n_perm)
    _, pa, pb, Dp = cp.observed()
    x = npm[np.isfinite(npm)]
    desc.append(dict(descriptive="BRCA D_PAM50 = pairing(luminal A+B) - pairing(basal + HER2)", value=Dp, p=float((1 + np.sum(x >= Dp)) / (len(x) + 1)),
                     null_sd=float(x.std(ddof=1)), n_pairs=int(sum(cp.n["BRCA"])), note=f"luminal {pa:.4f} (n {cp.n['BRCA'][0]}), basal+HER2 {pb:.4f} (n {cp.n['BRCA'][1]})"))
    pd.DataFrame(desc).to_csv(out / "descriptives.csv", index=False)
    print(t.round(4).T.to_string(), flush=True)
    print(pd.DataFrame(desc).round(4).to_string(), flush=True)
    log("tests stage done")


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "tests"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-sets", type=int, default=200)
    args = ap.parse_args(argv)
    cli.setup(args)
    out = HERE / "output" / ("" if args.n_perm == 10_000 and args.n_sets == 200 else "debug")
    out.mkdir(parents=True, exist_ok=True)
    (stage_power if args.stage == "power" else stage_tests)(args, out)


if __name__ == "__main__":
    main()
