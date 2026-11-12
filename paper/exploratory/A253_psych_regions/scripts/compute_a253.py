"""A253: are the clock relations between brain regions of one person changed in major depression or psychosis?

Registered in paper/docs/registered_experiments.md (A253). Per cohort, region pair and group: pairing r (|complex correlation|
of the two regions' clock-plane states over donors, states centred within region x group) and offset (arg of the
cross-moment). R_g = mean r over pairs (weights donors); D = R_control - R_case; A = mean |offset_case - offset_control|
(weights case donors); families pooled over cohorts with weights = case donors. Null: diagnosis labels permuted within cohort.

Stages
  power  design, P253.0 (controls only), null SD and MDE of D and A per family -> output/{design,power,precondition}.csv
  score  (after the power note is dated) everything registered -> output/score_{tests,pairs,descriptive}.csv
  `--placebo`: diagnosis labels (cases of the family and controls) shuffled once within cohort (code test; output/placebo/).

    python compute_a253.py --stage power|score [--smoke] [--placebo]

Choices not fixed by the registration
- Loaders: Ensembl ids mapped with the gene table of GSE231509, arrays with GPL570 (as A252: the row with the highest median
  per symbol); y = ln(count + 1) for counts, log2 x ln 2 for GSE53987, the deposited log FPKM x ln 2 for GSE102556.
- GSE53987 is deposited on the linear scale (median above 30) and is log2-transformed.
- A region pair enters a group's statistics with >= 8 donors having both regions; a cohort's D and A use the pairs present in
  both groups. Centring uses all donors of the group with that region.
- P253.0 null: control donors permuted independently within each region (10,000). Companion sets are drawn per cohort from the
  mean log expression of its control samples; set j is pooled across cohorts.
- Family p (one-sided): share of label permutations with D (or A) >= observed.
- The placebo shuffles labels with one fixed seed; psychosis and depression families are shuffled separately.
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
SCOUT = PAPER / "exploratory/scouts_2026-11-12/output"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A252_elapsed_time/scripts"))
sys.path.insert(0, str(PAPER / "exploratory/scouts_2026-11-12/scripts"))
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA  # noqa: E402
import compute_a252 as E  # noqa: E402  (collapse, ens_map, annot, states)
from gsm_meta import read_meta  # noqa: E402

OUT = HERE.parent / "output"
R = DATA / "psych_multiregion"
SEED, N_PERM, N_SETS, N_BOOT, MDE_K, MIN_N = 253, 10_000, 200, 2000, 2.49, 8
WF = cs.gtex()
HOURS = 24 / (2 * np.pi)
GROUP = {"Control": "control", "CTRL": "control", "control": "control", "match control": "control", "Major Depression": "mdd", "MDD": "mdd",
         "major depressive disorder": "mdd", "Schizophrenia": "scz", "schizophrenia": "scz", "psychosis_schizophrenia": "scz",
         "Bipolar Disorder": "bip", "bipolar disorder": "bip", "psychosis_bipolar": "bip"}
FAM = {"depression": (("mdd",), ["GSE80655", "GSE102556", "GSE53987"]), "psychosis": (("scz", "bip"), ["GSE80655", "GSE202537", "GSE53987"])}
COHORTS = ["GSE80655", "GSE102556", "GSE53987", "GSE202537"]
say = E.say


def load(acc):
    """(y: samples x genes log expression, meta: donor, region, group, tod)."""
    cache = HERE.parent / "output" / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    sh = pd.read_csv(SCOUT / f"psych_{acc}_samples.csv", dtype={"donor_id": str})
    sh["g"] = sh.group.map(GROUP)
    if acc == "GSE102556":
        mc = pd.read_csv(SCOUT / "psych_GSE102556_matrix_columns.csv", dtype=str)
        dg = sh.drop_duplicates("donor_id").set_index("donor_id").g
        meta = pd.DataFrame({"donor": mc.donor.to_numpy(), "region": mc.rcode.to_numpy(), "g": mc.donor.map(dg).to_numpy()}, index=mc.donor + "." + mc.rcode)
    else:
        key = {"GSE80655": "library_id", "GSE202537": "library_id", "GSE53987": "sample_id"}[acc]
        meta = pd.DataFrame({"donor": sh.donor_id.to_numpy(), "region": sh.region.to_numpy(), "g": sh.g.to_numpy(),
                             "tod": sh.tod_value.to_numpy() if "tod_value" in sh else np.nan}, index=sh[key].astype(str))
    f = cache / f"{acc}.parquet"
    if f.exists():
        y = pd.read_parquet(f)
    else:
        if acc == "GSE53987":
            e = read_meta(str(R / acc / "GSE53987_series_matrix.txt.gz"))[1]
            if float(np.nanmedian(e.to_numpy())) > 30:   # deposited on the linear scale
                e = np.log2(e.clip(lower=1.0))
            e = E.collapse(e, E.annot("GPL570").set_index("ID")["Gene symbol"]) * np.log(2.0)
        else:
            name = {"GSE80655": "GSE80655_GeneExpressionData_Updated_3-26-2018.txt.gz", "GSE202537": "GSE202537_NAc_Cau_Put_counts_corrected.csv.gz",
                    "GSE102556": "GSE102556_HumanMDD_fpkmtab.txt.gz"}[acc]
            e = pd.read_csv(R / acc / name, sep="," if name.endswith("csv.gz") else "\t", index_col=0)
            e = e.select_dtypes("number")
            e.index = e.index.astype(str).str.split(".").str[0]
            e = E.collapse(e[~e.index.duplicated()], E.ens_map())
            e = e * np.log(2.0) if acc == "GSE102556" else np.log(e + 1.0)
        y = e.T.astype(np.float32)
        y.index = y.index.astype(str)
        y.to_parquet(f)
    meta = meta[meta.index.isin(y.index) & meta.g.notna()]
    meta = meta[~meta.index.duplicated()]
    return y.loc[meta.index], meta


class Cohort:
    def __init__(self, acc):
        self.acc = acc
        self.y, self.meta = load(acc)
        self.regions = sorted(self.meta.region.unique())
        self.donors = sorted(self.meta.donor.unique())
        self.g = self.meta.drop_duplicates("donor").set_index("donor").g.reindex(self.donors).to_numpy()
        self.pairs = list(combinations(range(len(self.regions)), 2))
        self.row = pd.Index(self.donors).get_indexer(self.meta.donor)
        self.col = pd.Index(self.regions).get_indexer(self.meta.region)

    def Z(self, wf):
        z = E.states(self.y, wf)
        M = np.full((len(self.donors), len(self.regions)), np.nan + 0j)
        M[self.row, self.col] = z
        return M

    def pseudo(self, kind, n, tag):
        lm = self.y[(self.meta.g == "control").to_numpy()].mean(0).astype(float)
        g = [x for x in WF.genes if x in self.y.columns]
        return P.pseudo_waveforms(lm[~lm.index.duplicated()], n, np.random.default_rng([SEED, tag, sum(map(ord, self.acc))]), WF.subset(g), kind)


def group_stats(M, mask, pairs, perm_regions=None):
    """pair -> (r, cross-moment, n) for the donors in `mask`; states centred within region over those donors."""
    X = M[mask]
    X = X - np.nanmean(X, axis=0)
    if perm_regions is not None:
        X = np.column_stack([X[perm_regions[j], j] for j in range(X.shape[1])])
    out = {}
    for a, b in pairs:
        ok = ~(np.isnan(X[:, a]) | np.isnan(X[:, b]))
        if ok.sum() < MIN_N:
            continue
        c = np.sum(X[ok, a] * np.conj(X[ok, b]))
        out[(a, b)] = (abs(c) / np.sqrt(np.sum(abs(X[ok, a]) ** 2) * np.sum(abs(X[ok, b]) ** 2)), c, int(ok.sum()))
    return out


def R_of(st):
    w = np.array([v[2] for v in st.values()], float)
    return float(np.sum(w * [v[0] for v in st.values()]) / w.sum())


def contrast(M, lab, pairs):
    """(D, A, case donors) of one cohort for labels lab in {'control', 'case', other}."""
    sc, sk = group_stats(M, lab == "control", pairs), group_stats(M, lab == "case", pairs)
    common = [p for p in sc if p in sk]
    D = R_of({p: sc[p] for p in common}) - R_of({p: sk[p] for p in common})
    w = np.array([sk[p][2] for p in common], float)
    d = np.array([abs((np.angle(sk[p][1] * np.conj(sc[p][1])) * HOURS + 12) % 24 - 12) for p in common])
    return D, float(np.sum(w * d) / w.sum()), int((lab == "case").sum())


def family(Ms, labs, coh):
    v = [contrast(Ms[a], labs[a], coh[a].pairs) for a in labs]
    w = np.array([x[2] for x in v], float)
    return float(np.sum(w * [x[0] for x in v]) / w.sum()), float(np.sum(w * [x[1] for x in v]) / w.sum())


def labels(c, cases):
    return np.where(c.g == "control", "control", np.where(np.isin(c.g, cases), "case", "other"))


def perm_labels(lab, rng):
    out = lab.copy()
    i = np.flatnonzero(lab != "other")
    out[i] = lab[rng.permutation(i)]
    return out


def p_up(nul, obs):
    return float((1 + np.sum(np.asarray(nul) >= obs)) / (1 + len(nul)))


def holm2(a, b):
    lo, hi = sorted([a, b])
    return {a: min(1, 2 * lo) if a == lo else max(min(1, 2 * lo), hi), b: min(1, 2 * lo) if b == lo else max(min(1, 2 * lo), hi)}


def run(args):
    n_perm, n_sets, n_boot = (300, 5, 100) if args.smoke else (N_PERM, N_SETS, N_BOOT)
    coh = {a: Cohort(a) for a in COHORTS}
    Z = {a: coh[a].Z(WF) for a in COHORTS}
    design = []
    for a, c in coh.items():
        for g in pd.unique(c.g):
            st = group_stats(Z[a], c.g == g, c.pairs)
            design.append(dict(cohort=a, group=g, donors=int((c.g == g).sum()), regions=len(c.regions), pairs=len(st),
                               donors_per_pair_min=min((v[2] for v in st.values()), default=0), clock_genes=len([x for x in WF.genes if x in c.y.columns])))
    pd.DataFrame(design).to_csv(OUT / "design.csv", index=False)
    lab = {f: {a: labels(coh[a], cases) for a in cs_} for f, (cases, cs_) in FAM.items()}
    if args.placebo:
        rng = np.random.default_rng([SEED, 9])
        lab = {f: {a: perm_labels(v, rng) for a, v in d.items()} for f, d in lab.items()}

    # P253.0: controls
    ctrl = {a: coh[a].g == "control" for a in COHORTS}
    def Rc(Ms, perm=None):
        st = [group_stats(Ms[a], ctrl[a], coh[a].pairs, None if perm is None else perm[a]) for a in COHORTS]
        w = np.array([ctrl[a].sum() for a in COHORTS], float)
        return float(np.sum(w * [R_of(s) for s in st]) / w.sum())
    rng = np.random.default_rng([SEED, 1])
    obs0 = Rc(Z)
    nul0 = np.array([Rc(Z, {a: [rng.permutation(int(ctrl[a].sum())) for _ in coh[a].regions] for a in COHORTS}) for _ in range(n_perm)])
    sets = {kind: {a: coh[a].pseudo(kind, n_sets, tag) for a in COHORTS} for kind, tag in (("matched", 3), ("responsive", 4))}
    Zs = {kind: [{a: coh[a].Z(sets[kind][a][j]) for a in COHORTS} for j in range(n_sets)] for kind in sets}
    pre = dict(R_control=obs0, null_mean=float(nul0.mean()), null_sd=float(nul0.std()), p=p_up(nul0, obs0))
    for kind in sets:
        v = np.array([Rc(m) for m in Zs[kind]])
        pre.update({f"{kind}_median": float(np.median(v)), f"{kind}_q95": float(np.quantile(v, .95)), f"p_{kind}": p_up(v, obs0)})
    pre["passes"] = bool(pre["p"] < 0.05 and obs0 > pre["matched_q95"] and obs0 > pre["responsive_q95"])
    pd.DataFrame([pre]).to_csv(OUT / "precondition.csv", index=False)

    nul = {}
    rng = np.random.default_rng([SEED, 2])
    for f in FAM:
        nul[f] = np.array([family(Z, {a: perm_labels(v, rng) for a, v in lab[f].items()}, coh) for _ in range(n_perm)])
    power = pd.DataFrame([dict(family=f, cohorts=len(lab[f]), cases=int(sum((v == "case").sum() for v in lab[f].values())),
                               controls=int(sum((v == "control").sum() for v in lab[f].values())), D_null_mean=float(nul[f][:, 0].mean()),
                               D_null_sd=float(nul[f][:, 0].std()), D_mde=MDE_K * float(nul[f][:, 0].std()), A_null_mean=float(nul[f][:, 1].mean()),
                               A_null_sd=float(nul[f][:, 1].std()), A_mde_h=MDE_K * float(nul[f][:, 1].std())) for f in FAM])
    power["scored"] = bool(pre["passes"] and power.D_mde.iloc[0] <= 0.20)
    power.to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(design).to_string(), pd.Series(pre).round(4).to_string(), power.round(4).to_string(), sep="\n", flush=True)
    if args.stage == "power":
        return

    tests, desc, prow = [], [], []
    rng = np.random.default_rng([SEED, 5])
    for f in FAM:
        D, A = family(Z, lab[f], coh)
        boot = []
        for _ in range(n_boot):
            Mb, lb = {}, {}
            for a, v in lab[f].items():
                i = np.concatenate([rng.choice(np.flatnonzero(v == g), int((v == g).sum())) for g in ("control", "case")])
                Mb[a], lb[a] = Z[a][i], v[i]
            boot.append(family(Mb, lb, coh)[0])
        for name, obs, col in (("offsets (A, h)", A, 1), ("pairing (D)", D, 0)):
            row = dict(family=f, test=name, value=obs, null_mean=float(nul[f][:, col].mean()), null_sd=float(nul[f][:, col].std()),
                       p=p_up(nul[f][:, col], obs), mde=MDE_K * float(nul[f][:, col].std()))
            for kind in sets:
                v = np.array([family(m, lab[f], coh)[col] for m in Zs[kind]])
                row.update({f"{kind}_median": float(np.median(v)), f"{kind}_q95": float(np.quantile(v, .95)), f"p_{kind}": p_up(v, obs)})
            row["clock_specific"] = bool(obs > row["matched_q95"] and obs > row["responsive_q95"])
            if col == 0:
                row.update(D_lo90=float(np.quantile(boot, .05)), D_hi90=float(np.quantile(boot, .95)))
            tests.append(row)
        for a, v in lab[f].items():
            d, aa, n = contrast(Z[a], v, coh[a].pairs)
            desc.append(dict(item=f"{f}: per cohort", cohort=a, D=d, A_h=aa, cases=n, controls=int((v == "control").sum())))
    T = pd.DataFrame(tests)
    for f in FAM:
        m = T.family == f
        h = holm2(*T.loc[m, "p"])
        T.loc[m, "p_holm"] = T.loc[m, "p"].map(h)
    for a, c in coh.items():
        for g in pd.unique(c.g):
            m = c.g == g
            st = group_stats(Z[a], m, c.pairs)
            for (i, j), (r, cm, n) in st.items():
                prow.append(dict(cohort=a, group=g, region_a=c.regions[i], region_b=c.regions[j], donors=n, r=float(r), offset_h=float(np.angle(cm) * HOURS)))
            X = Z[a][m] - np.nanmean(Z[a][m], axis=0)
            desc.append(dict(item="clock-plane variance (mean over regions)", cohort=a, group=g, value=float(np.nanmean(abs(X) ** 2)), R=R_of(st) if st else np.nan, donors=int(m.sum())))
        if a in ("GSE80655", "GSE53987", "GSE202537"):
            for dx in ("scz", "bip"):
                v = np.where(c.g == "control", "control", np.where(c.g == dx, "case", "other"))
                if (v == "case").sum() >= MIN_N:
                    d, aa, n = contrast(Z[a], v, c.pairs)
                    desc.append(dict(item=f"{dx} alone", cohort=a, D=d, A_h=aa, cases=n))
    c = coh["GSE202537"]
    tod = pd.Series(c.meta.tod.to_numpy(float), index=c.meta.donor).groupby(level=0).first().reindex(c.donors).to_numpy()
    for g in pd.unique(c.g):
        for j, reg in enumerate(c.regions):
            m = (c.g == g) & ~np.isnan(Z["GSE202537"][:, j]) & ~np.isnan(tod)
            z = Z["GSE202537"][m, j]
            a = np.mean((z - z.mean()) * np.exp(-1j * tod[m] / HOURS))
            desc.append(dict(item="GSE202537 cosinor of state on time of death", group=g, region=reg, donors=int(m.sum()), value=float(abs(a)), acrophase_shift_h=float(np.angle(a) * HOURS)))
    T.to_csv(OUT / "score_tests.csv", index=False)
    pd.DataFrame(prow).to_csv(OUT / "score_pairs.csv", index=False)
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    print(T.round(4).T.to_string(), pd.DataFrame(desc).round(3).to_string(), sep="\n", flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--placebo", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args(argv)
    cli.setup(args)
    global OUT
    OUT = OUT / ("placebo" if args.placebo else "smoke" if args.smoke else "")
    OUT.mkdir(parents=True, exist_ok=True)
    run(args)


if __name__ == "__main__":
    main()
