"""A229: do tumour endothelial adhesion genes and the tumour CD8 T share follow the tumour endothelial clock state?

Registered in paper/docs/registered_experiments.md (A229). The 34 A227 cohorts and tables (`compute_a227`), A227 / A149 states and
pseudo-clock construction, readouts y centred within cohort x chemistry (the strata of the states).

State z: complex clock state of the endothelial pseudobulk of each patient (`clocksync.paired.projector / states`, GTEx waveform,
counts=True), projector fitted on the cohort's endothelial rows of that side only, centred within stratum over the patients of the test.
Statistic: r = |sum y conj(z)| / sqrt(sum y^2 sum |z|^2) per cohort (y real, centred), pooled with patient-count weights.
Null: y permuted across patients within cohort. Clock specificity: z from expression-matched pseudo-clock waveforms (A227 `pseudo_one`
construction, reference = the endothelial genome pseudobulks), same y; p_clock = (1 + #{pseudo >= obs}) / (1 + n_pseudo).

Tests (y): H229.1 ICAM1, H229.3 VCAM1, H229.4 SELE (log2(CPM + 1) of the genome-wide endothelial pseudobulk); H229.2 logit CD8 share
log((CD8T + .5) / (CD8neg_T + .5)) of the T cells (A228 `subset_fractions.csv`, samples of the side, summed per patient, >= 50 T/NK cells);
H229.5 log of the T/NK share of all cells of the sample (n_cells of the T_NK pseudobulk over the sum over the compartment pseudobulks of the
table, so compartments below the A149 rules are not counted). Primaries H229.1-2, secondaries H229.3-5 (Holm separately).
Side `tumour` is the registered test; `adjacent_normal` runs the same tests as descriptive.

  power  eligibility and the within-cohort permutation null from SHUFFLED y only (the identity pairing is never evaluated; no observed
         statistic is computed, stored or printed); MDE = 2.49 x pooled null SD; scored = >= 40 patients and MDE <= 0.25.
                                                                    -> output/{eligibility,power}.csv
  tests  (IMPLEMENTED, NOT RUN before the power note is committed) observed pooled r, permutation p, pseudo q95, Holm, pass; descriptive:
         ICAM1 percentile among 200 expression-matched endothelial genes as y, adjacent-normal side, cancer types with >= 12 patients, implied
         peak hour.                          -> output/{tests,pooled,per_cohort,pseudo,matched_genes}.csv

Implied peak: y is high where the state points to angle theta = arg(sum y z) = -arg(sum y conj z); hours = (theta mod 2 pi) * 24 / (2 pi)
(`clocksync.stats.HOURS`, the convention of `Dataset.host` `phase_h`: angle of the complex clock state on the GTEx waveform axis).

    python compute_a229.py --stage power|tests [--jobs J] [--n-perm 10000] [--n-pseudo 200] [--force] [--smoke --out DIR]
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
for p in ("A227_tme_compartment_pairing", "A149_paired_tumour_normal", "A151_fibroblast_pairing", "A226_stroma_pairing_confirm"):
    sys.path.insert(0, str(PAPER / "exploratory" / p / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
import compute_a227 as A  # noqa: E402
from clocksync import stats as S  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
FRACTIONS = PAPER / "exploratory" / "A228_immune_body_time" / "output" / "subset_fractions.csv"
SEED = 229
SIDES = ["tumour", "adjacent_normal"]
# test -> (hypothesis, role, kind of readout, gene)
TESTS = {"ICAM1": ("H229.1", "primary", "gene", "ICAM1"), "CD8_share": ("H229.2", "primary", "cd8", None),
         "VCAM1": ("H229.3", "secondary", "gene", "VCAM1"), "SELE": ("H229.4", "secondary", "gene", "SELE"),
         "TNK_share": ("H229.5", "secondary", "tnk", None)}
KINDS = list(TESTS)
GENES = sorted({v[3] for v in TESTS.values() if v[3]})
MIN_PATIENTS, MDE_K, MAX_MDE, MIN_CANCER, MIN_N, MIN_TNK = A.MIN_PATIENTS, A.MDE_K, A.MAX_MDE, A.MIN_CANCER, L.MIN_PAIRED, 50
N_MATCH, N_POOL = 200, 400


def seed_of(kind, study, side):
    return [SEED, KINDS.index(kind), SIDES.index(side), sum(map(ord, study))]


# ---- per-cohort data ---------------------------------------------------------------------------------------------

def _comp(side):
    return f"{side}:endothelial"


def cohort_data(study, side, n_pseudo=0, full=False):
    """Endothelial states and readout sources of one cohort and side. None if < MIN_N endothelial patients.
    Keys: pts, stratum, z (uncentred states, real waveform), lcpm (patients x genes log2(CPM+1); the test genes only unless `full`),
    Zp (n_pseudo x patients pseudo states), cells (all-compartment cell count), tnk (T_NK cells), n_CD8T, n_CD8neg_T, n_T_NK."""
    from clocksync.config import limit_threads
    from clocksync.paired import projector, pseudo_waveforms, states
    limit_threads(1)
    pb = A.table(study)
    comp = _comp(side)
    rows = L._rows(pb, comp)
    if len(rows) < MIN_N:
        return dict(study=study, side=side, pts=[], reason=f"{len(rows)} endothelial patients")
    G = A.genome_pair(study, (comp,), list(rows.index))
    G = G[G.compartment == comp].drop_duplicates("donor").set_index("donor")
    pts = sorted(set(rows.index) & set(G.index))
    if len(pts) < MIN_N:
        return dict(study=study, side=side, pts=[], reason=f"{len(pts)} endothelial patients in the genome table")
    rows, G = rows.loc[pts], G.loc[pts]
    z = S.complex_state(states(rows, projector(rows, L.q.WF, True), L.q.WF, depth=rows.depth.to_numpy(float))[0])
    genes = [c for c in G.columns if c not in ("compartment", "stratum", "depth")]
    X = G[genes].to_numpy(float)
    lc = np.log2(X / G.depth.to_numpy(float)[:, None] * 1e6 + 1)
    keep = np.ones(len(genes), bool) if full else np.isin(genes, [L.canonical(g) for g in GENES])
    lcpm = pd.DataFrame(lc[:, keep].astype(np.float32), index=pts, columns=np.asarray(genes)[keep])
    Zp = np.zeros((0, len(pts)), complex)
    if n_pseudo:
        from clocksync.paired import pseudo_waveforms
        frac = G[genes].div(G.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(G[genes].mean() >= 1) | frac.columns.isin(list(L.q.WF.genes))]
        rng = np.random.default_rng(seed_of("ICAM1", study, side) + [44])
        Zp = np.array([S.complex_state(states(G, projector(G, w, True), w, depth=G.depth.to_numpy(float))[0])
                       for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF)])
    pbk = pb[pb.compartment.str.startswith(side + ":")].drop_duplicates(["donor", "compartment"])
    cells = pbk.groupby("donor").n_cells.sum().reindex(pts)
    tnk = pbk[pbk.compartment == f"{side}:T_NK"].set_index("donor").n_cells.reindex(pts)
    d = dict(study=study, side=side, pts=pts, stratum=rows.stratum.astype(str).to_numpy(), z=z, lcpm=lcpm, Zp=Zp, cells=cells, tnk=tnk)
    fr = pd.read_csv(FRACTIONS, dtype={"patient": str, "sample": str})
    fr = fr[(fr.study == study) & (fr.kind == side)].groupby("patient")[["n_T_NK", "n_CD8T", "n_CD8neg_T"]].sum()
    for c in fr:
        d[c] = fr[c].reindex(pts)
    d["has_fractions"] = len(fr) > 0
    return d


def readout(d, kind):
    """(y over d['pts'] with NaN where undefined, reason if the whole cohort drops out)."""
    _, _, how, gene = TESTS[kind]
    n = len(d["pts"])
    if how == "gene":
        g = L.canonical(gene)
        return (d["lcpm"][g].to_numpy(float), "") if g in d["lcpm"].columns else (np.full(n, np.nan), f"{gene} absent")
    if how == "cd8":
        if not d["has_fractions"]:
            return np.full(n, np.nan), "no subset fractions (A228)"
        a, b, t = d["n_CD8T"].to_numpy(float), d["n_CD8neg_T"].to_numpy(float), d["n_T_NK"].to_numpy(float)
        return np.where(t >= MIN_TNK, np.log((a + .5) / (b + .5)), np.nan), ""
    t, c = d["tnk"].to_numpy(float), d["cells"].to_numpy(float)
    return np.where(t > 0, np.log(t / c), np.nan), ""


def _centre(v, stratum):
    return v - pd.Series(v).groupby(stratum).transform("mean").to_numpy()


def make_unit(d, y=None, kind=None):
    """Centred (within stratum, over the patients with a defined y) y and z; None if < MIN_N patients or a degenerate spread."""
    y, reason = readout(d, kind) if y is None else (y, "")
    m = np.isfinite(y)
    if m.sum() < MIN_N:
        return None, reason or f"{int(m.sum())} patients with the readout"
    st = d["stratum"][m]
    yc, zc = _centre(y[m], st), _centre(d["z"][m], st)
    if not (np.sum(yc ** 2) > 0 and np.sum(abs(zc) ** 2) > 0):
        return None, "no spread in y or z"
    return dict(n=int(m.sum()), mask=m, stratum=st, y=yc, z=zc, den=float(np.sqrt(np.sum(yc ** 2) * np.sum(abs(zc) ** 2)))), ""


def cross(y, z):
    """sum y conj(z) of centred vectors (complex)."""
    return np.sum(y * np.conj(z))


def stat(u, y=None, z=None):
    y, z = u["y"] if y is None else y, u["z"] if z is None else z
    return float(abs(cross(y, z)) / np.sqrt(np.sum(y ** 2) * np.sum(abs(z) ** 2)))


def perm_null(u, n_perm, rng):
    """Null: y shuffled across patients (n_perm uniform permutations). The identity pairing is not evaluated."""
    perms = np.argsort(rng.random((n_perm, u["n"])), axis=1)
    return np.abs((u["y"][perms] * np.conj(u["z"])[None, :]).sum(1)) / u["den"]


def pseudo_stats(d, u):
    """r of each pseudo-clock set with the same y (centred within the same strata, same patients)."""
    return np.array([stat(u, z=_centre(z[u["mask"]], u["stratum"])) for z in d["Zp"]])


def peak_h(c):
    return float((-np.angle(c)) % (2 * np.pi) * S.HOURS)


# ---- power -----------------------------------------------------------------------------------------------------

def power_one(study, side, n_perm):
    d = cohort_data(study, side)
    elig, nulls = [], {}
    for kind in KINDS:
        base = dict(study=study, cancer=A.cancer(study), side=side, kind=kind, hypothesis=TESTS[kind][0],
                    patients_endothelial=len(d["pts"]))
        if not d["pts"]:
            elig.append(dict(base, patients=0, eligible=False, reason=d["reason"]))
            continue
        u, why = make_unit(d, kind=kind)
        elig.append(dict(base, patients=u["n"] if u else 0, eligible=u is not None, reason=why))
        if u:
            nulls[kind] = (u["n"], perm_null(u, n_perm, np.random.default_rng(seed_of(kind, study, side))))
    return elig, (study, side, nulls)


def stage_power(args):
    cohorts = A.cohort_list(args.smoke)
    res = Parallel(n_jobs=args.jobs)(delayed(power_one)(s, sd, args.n_perm) for s in cohorts for sd in SIDES)
    elig = pd.DataFrame([r for e, _ in res for r in e])
    elig.to_csv(args.out / "eligibility.csv", index=False)
    nl = {(s, sd): n for _, (s, sd, n) in res}
    rows = []
    for side in SIDES:
        for kind in KINDS:
            cs = [s for s in cohorts if kind in nl[(s, side)]]
            ns = [nl[(s, side)][kind][0] for s in cs]
            nul = A.pool_null(ns, [nl[(s, side)][kind][1] for s in cs]) if cs else np.array([np.nan])
            sd_, n = float(nul.std()), int(sum(ns))
            rows.append(dict(side=side, kind=kind, hypothesis=TESTS[kind][0], role=TESTS[kind][1] if side == "tumour" else "descriptive",
                             patients=n, cohorts=len(cs), cohort_list=";".join(cs), null_mean=float(nul.mean()), null_sd=sd_,
                             mde=MDE_K * sd_, n_perm=args.n_perm, scored=bool(n >= MIN_PATIENTS and MDE_K * sd_ <= MAX_MDE)))
    P = pd.DataFrame(rows)
    P.to_csv(args.out / "power.csv", index=False)
    print(P.drop(columns="cohort_list").round(4).to_string(), flush=True)


# ---- tests (registered statistic; not run before the power note) ----------------------------------------------

def stage_tests(args):
    pw = pd.read_csv(args.out / "power.csv")
    pw = pw[pw.side == "tumour"].set_index("kind")
    cohorts = A.cohort_list(args.smoke)
    kinds = [k for k in KINDS if bool(pw.loc[k, "scored"]) or args.force]
    for k in KINDS:
        if k not in kinds:
            print(f"{k}: not scored (power note); use --force to run anyway", flush=True)
    data = Parallel(n_jobs=args.jobs)(delayed(cohort_data)(s, sd, args.n_pseudo, True) for s in cohorts for sd in SIDES)
    data = {(d["study"], d["side"]): d for d in data if d["pts"]}
    tests, pooled, per, pseudo_rows, matched, icam_pct = [], [], [], [], [], np.nan

    def pool(items, kind, side):
        """items: [(study, d, u)]. Pooled observed r, null, pseudo, complex cross-correlation, weights."""
        w = np.array([u["n"] for _, _, u in items], float)
        obs = float(np.sum(w * [stat(u) for _, _, u in items]) / w.sum())
        c = np.sum(w * [cross(u["y"], u["z"]) / u["den"] for _, _, u in items]) / w.sum()
        nul = A.pool_null(w, [perm_null(u, args.n_perm, np.random.default_rng(seed_of(kind, s, side) + [1])) for s, _, u in items])
        pse = A.pooled_pseudo(np.array([pseudo_stats(d, u) for _, d, u in items]), w)
        return w, obs, nul, pse, c

    for side in SIDES:
        for kind in kinds:
            items = []
            for s in cohorts:
                if (s, side) in data:
                    u, _ = make_unit(data[(s, side)], kind=kind)
                    if u:
                        items.append((s, data[(s, side)], u))
            if not items:
                continue
            w, obs, nul, pse, c = pool(items, kind, side)
            p1, p2 = (1 + np.sum(nul >= obs)) / (1 + len(nul)), (1 + np.sum(pse >= obs)) / (1 + len(pse))
            tests.append(dict(side=side, kind=kind, hypothesis=TESTS[kind][0], role=TESTS[kind][1] if side == "tumour" else "descriptive",
                              patients=int(w.sum()), cohorts=len(items), statistic=obs, null_mean=float(nul.mean()), null_sd=float(nul.std()),
                              p_pairing=float(p1), pseudo_median=float(np.median(pse)), pseudo_q95=float(np.quantile(pse, .95)),
                              n_sets=len(pse), p_clock=float(p2), clock_minus_pseudo_median=obs - float(np.median(pse)),
                              arg_sum_y_conj_z_h=float(np.angle(c) % (2 * np.pi) * S.HOURS), implied_peak_h=peak_h(c)))
            pooled.append(tests[-1].copy())
            pseudo_rows += [dict(side=side, kind=kind, set=i, value=v) for i, v in enumerate(pse)]
            for s, d, u in items:
                per.append(dict(side=side, kind=kind, cohort=s, cancer=A.cancer(s), n=u["n"], statistic=stat(u)))
            for ca in sorted({A.cancer(s) for s, _, _ in items}):          # descriptive: cancer types with >= 12 patients
                sc = [x for x in items if A.cancer(x[0]) == ca]
                if sum(u["n"] for _, _, u in sc) >= MIN_CANCER:
                    wc, oc, _, pc, cc = pool(sc, kind, side)
                    pooled.append(dict(side=side, kind=kind, cancer=ca, patients=int(wc.sum()), cohorts=len(sc), statistic=oc,
                                       pseudo_median=float(np.median(pc)), clock_minus_pseudo_median=oc - float(np.median(pc)),
                                       implied_peak_h=peak_h(cc)))
            if side == "tumour" and kind == "ICAM1":
                matched, icam_pct = icam_percentile(items, tests[-1], args)
    Tt = pd.DataFrame(tests)
    Tt["p_pairing_holm"] = np.nan
    for role in ("primary", "secondary"):
        m = (Tt.side == "tumour") & (Tt.role == role)
        if m.any():
            Tt.loc[m, "p_pairing_holm"] = A.holm(Tt.loc[m, "p_pairing"])
    Tt["pass"] = np.where(Tt.side == "tumour", (Tt.p_pairing_holm < 0.05) & (Tt.statistic > Tt.pseudo_q95), np.nan)
    Tt["icam1_percentile_among_matched"] = np.where((Tt.side == "tumour") & (Tt.kind == "ICAM1"), icam_pct, np.nan)
    if matched:
        pd.DataFrame(matched).to_csv(args.out / "matched_genes.csv", index=False)
    Tt.to_csv(args.out / "tests.csv", index=False)
    pd.DataFrame(pooled).to_csv(args.out / "pooled.csv", index=False)
    pd.DataFrame(per).to_csv(args.out / "per_cohort.csv", index=False)
    pd.DataFrame(pseudo_rows).to_csv(args.out / "pseudo.csv", index=False)
    print(Tt.round(4).to_string(), flush=True)


def icam_percentile(items, test_row, args):  # noqa: ARG001
    """Percentile of ICAM1's r among N_MATCH expression-matched endothelial genes used as y on the same patients and states.
    Genes present in all cohorts of the ICAM1 test, excluding the test genes and the waveform genes; the N_MATCH are drawn (seeded) from
    the N_POOL genes closest to ICAM1 in patient-weighted mean log2(CPM + 1); r pooled as the test."""
    names = set.intersection(*[set(d["lcpm"].columns) for _, d, _ in items])
    excl = {L.canonical(g) for g in GENES} | set(L.q.WF.genes)
    w = np.array([u["n"] for _, _, u in items], float)
    mean = pd.Series({g: float(np.sum(w * [d["lcpm"][g].to_numpy(float)[u["mask"]].mean() for _, d, u in items]) / w.sum())
                      for g in sorted(names)})
    ic = mean[L.canonical("ICAM1")]
    cand = (mean.drop([g for g in excl if g in mean.index]) - ic).abs().sort_values().index[:N_POOL]
    pick = np.random.default_rng([SEED, 7]).choice(cand, size=min(N_MATCH, len(cand)), replace=False)
    out = []
    for g in pick:
        rs, ws = [], []
        for _, d, u in items:
            yc = _centre(d["lcpm"][g].to_numpy(float)[u["mask"]], u["stratum"])
            if np.sum(yc ** 2) > 0:
                rs.append(stat(u, y=yc))
                ws.append(u["n"])
        out.append(dict(gene=g, mean_log2cpm=float(mean[g]), statistic=float(np.average(rs, weights=ws)) if rs else np.nan))
    r = np.array([o["statistic"] for o in out])
    pct = float(np.mean(r[np.isfinite(r)] <= test_row["statistic"]) * 100)
    for o in out:
        o["icam1_statistic"], o["icam1_mean_log2cpm"] = test_row["statistic"], float(ic)
    return out, pct


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "tests"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--force", action="store_true", help="tests: also run tests flagged not scored")
    ap.add_argument("--smoke", action="store_true", help="3 cohorts, 200 permutations, 5 pseudo sets")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    assert args.jobs <= 3, "at most 3 workers for A229"
    cli.setup(args)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.smoke:
        args.n_perm, args.n_pseudo = 200, 5
    {"power": stage_power, "tests": stage_tests}[args.stage](args)


if __name__ == "__main__":
    main()
