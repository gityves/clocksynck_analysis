"""A252: does the clock state of a living person's tissue advance by the elapsed time between two biopsies of one day?

Registered in paper/docs/registered_experiments.md (A252). Unit = person x visit (or leg). State z = L y (unweighted
projector of the GTEx waveform on the clock genes present; y = log2 x ln 2 for arrays, ln(count + 1) for counts).
For a unit and a later time point k: delta = z_k - z_0, predicted chord c_k = exp(i w (t0 + D_k)) - exp(i w t0);
rho = sum Re(delta conj(c)) / sum |c|^2 (pooled sums over datasets). Null: time-point labels permuted within unit.
Companions: rho from expression-matched and responsive pseudo-clock sets carrying the waveform's phases.

Stages
  power  units, clock genes, null SD and MDE of rho per family (permuted labels only)          -> output/power.csv
  score  (after the power note is dated) everything registered -> output/score_{tests,datasets,descriptive}.csv
  `--placebo`: time-point labels shuffled once within unit before anything else (code test; output/placebo/).

    python compute_a252.py --stage power|score [--smoke] [--placebo]

Choices not fixed by the registration
- Loaders follow `scouts_2026-11-12/scripts/build_sheets.py` (sample sheets from the scout). Gene symbols: arrays collapsed to
  the probe with the highest median; Ensembl ids mapped with the gene table of GSE231509; Entrez ids (GSE118280) with
  GPL570's annotation. Duplicated symbols: the first.
- Units need the first time point and >= 1 later one. GSE303983: elapsed 2.25 h for every person (arm not public).
- Pseudo-clock sets are drawn per dataset from its genome-wide mean log expression; set j is pooled across datasets.
- Family p: share of permutations with rho >= observed (10,000 random label permutations per unit).
- Rest vs exercise: D = mean over common later time points of |mean delta(exercise) - mean delta(rest)|; GSE252357: arm
  labels permuted over persons; GSE120862: legs swapped within person x training status.
- An array matrix with a median above 30 is taken as linear and log2-transformed (GSE4386).
- Descriptive datasets: modulus of the mean delta against the same modulus of the companion sets.
"""

from __future__ import annotations

import gzip
import io
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
SCOUT = PAPER / "exploratory/scouts_2026-11-12"
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(SCOUT / "scripts"))
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA  # noqa: E402
from gsm_meta import read_meta  # noqa: E402

OUT = HERE.parent / "output"
R = DATA / "repeat_biopsy"
SEED, N_PERM, N_SETS, MDE_K, MAX_MDE = 252, 10_000, 200, 2.49, 0.7
WF = cs.gtex()
W = 2 * np.pi / 24
# dataset -> (tissue, clock time of the first biopsy, counts?, role)
DS = {"GSE118280": ("adipose", 8 + 40 / 60, False, "primary"), "GSE231509": ("muscle", 9.5, True, "pooled"),
      "GSE252357": ("muscle", 9.0, True, "pooled"), "GSE142401": ("adipose", 9.0, False, "pooled"),
      "GSE151066": ("muscle", 9.0, True, "pooled"), "GSE273041": ("muscle", 9.0, True, "pooled"),
      "GSE303983": ("adipose", 9.0, True, "pooled"), "GSE120862": ("muscle", 9.0, True, "pooled"),
      "GSE4386": ("atrium", np.nan, False, "descriptive"), "GSE75432": ("muscle", np.nan, True, "descriptive"),
      "GSE151648": ("liver graft", np.nan, True, "descriptive")}
FAMILIES = {"H252.1 GSE118280 (exact times)": ["GSE118280"],
            "H252.2 seven datasets, assumed times": [d for d, v in DS.items() if v[3] == "pooled"],
            "H252.3 muscle": [d for d, v in DS.items() if v[3] == "pooled" and v[0] == "muscle"],
            "H252.3 adipose": [d for d, v in DS.items() if v[3] == "pooled" and v[0] == "adipose"]}


def say(*a):
    import time
    print(time.strftime("%F %T"), *a, flush=True)


def annot(gpl, col="Gene symbol"):
    lines = gzip.open(R / "platforms" / f"{gpl}.annot.gz", "rt", errors="replace").read().split("\n")
    i = next(k for k, l in enumerate(lines) if l.startswith("ID\t"))
    return pd.read_csv(io.StringIO("\n".join(lines[i:])), sep="\t", dtype=str, low_memory=False)


def collapse(expr, sym):
    """Rows of expr relabelled by `sym` (Series index -> symbol); the row with the highest median per symbol."""
    s = sym.reindex(expr.index)
    ok = s.notna().to_numpy() & (s.astype(str).str.strip() != "").to_numpy()
    e = expr[ok].copy()
    e.index = [canonical(str(x).split("///")[0].strip()) for x in s[ok]]
    order = np.argsort(-e.median(axis=1).to_numpy(), kind="stable")
    e = e.iloc[order]
    return e[~e.index.duplicated()]


def ens_map():
    t = pd.read_csv(R / "GSE231509/GSE231509_Read_counts.txt.gz", sep="\t", usecols=["gene_id", "gene_name"])
    return pd.Series(t.gene_name.to_numpy(), index=t.gene_id.str.split(".").str[0])


def matrix(acc):
    """genes (symbols) x samples (the sheet's matrix_column)."""
    if acc == "GSE252357":
        e = pd.read_csv(R / acc / "GSE252357_normalized_counts_protein_coding.csv.gz", index_col=0)
        return collapse(e, pd.Series(e.index, index=e.index))
    if acc in ("GSE118280", "GSE142401", "GSE4386"):
        e = read_meta(str(R / acc / f"{acc}_series_matrix.txt.gz"))[1]
        if acc == "GSE118280":
            a = annot("GPL570")
            m = {}
            for ids, syms in zip(a["Gene ID"].fillna(""), a["Gene symbol"].fillna("")):
                for i, s in zip(ids.split("///"), syms.split("///")):
                    m[i] = s
            e.index = e.index.astype(str).str.replace("_at", "")
            return collapse(e, pd.Series(m))
        a = annot("GPL10558" if acc == "GSE142401" else "GPL570")
        return collapse(e, a.set_index("ID")["Gene symbol"])
    if acc == "GSE231509":
        e = pd.read_csv(R / acc / "GSE231509_Read_counts.txt.gz", sep="\t", index_col=1).drop(columns=["gene_id", "gene_biotype"])
        return collapse(e, pd.Series(e.index, index=e.index).drop_duplicates())
    if acc == "GSE303983":
        e = pd.read_csv(R / acc / "GSE303983_genes.counts.csv.gz", index_col=0)
        e.index = e.index.str.split("|").str[-1].str.strip()
        e = e[~e.index.duplicated()]
        return collapse(e, pd.Series(e.index, index=e.index))
    if acc == "GSE120862":
        e = pd.read_csv(R / acc / "GSE120862_DESeq_normalized_counts.txt.gz", sep="\t", index_col=0)
        e = e[~e.index.duplicated()]
        return collapse(e, pd.Series(e.index, index=e.index))
    if acc == "GSE75432":
        cols = {}
        with tarfile.open(R / acc / "GSE75432_RAW.tar") as tf:
            for m in tf.getmembers():
                cols[m.name.split("_")[0]] = pd.read_csv(io.BytesIO(gzip.decompress(tf.extractfile(m).read())), sep="\t", header=None, index_col=0)[1]
        e = pd.DataFrame(cols)
    elif acc == "GSE151648":
        e = pd.read_csv(R / acc / "GSE151648_liver-iri-counts.txt.gz", sep="\t", index_col=0)
    elif acc == "GSE273041":
        e = pd.read_csv(R / acc / "GSE273041_raw_counts.txt.gz", sep="\t", index_col=0)
    elif acc == "GSE151066":
        e = pd.read_csv(R / acc / "GSE151066_rsem_genes_count.txt.gz", sep=r"\s+", quotechar='"')
    e = e[e.index.astype(str).str.startswith("ENSG")]
    e.index = e.index.astype(str).str.split(".").str[0]
    e = e[~e.index.duplicated()]
    return collapse(e, ens_map())


def load(acc):
    """(log expression samples x genes, sheet with unit / k / elapsed), cached."""
    cache = OUT.parents[0] / "output" / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    f = cache / f"{acc}.parquet"
    sh = pd.read_csv(SCOUT / "output" / f"repeat_{acc}_samples.csv")
    sh = sh[sh.matrix_column.notna()].copy()
    if f.exists():
        y = pd.read_parquet(f)
    else:
        e = matrix(acc).apply(pd.to_numeric, errors="coerce")
        e = e[[c for c in sh.matrix_column if c in e.columns]].dropna(how="any")
        if not DS[acc][2] and float(np.nanmedian(e.to_numpy())) > 30:   # array deposited on the linear scale (GSE4386)
            e = np.log2(e.clip(lower=1.0))
        y = (np.log(e + 1.0) if DS[acc][2] else e * np.log(2.0)).T.astype(np.float32)
        y.to_parquet(f)
    sh = sh[sh.matrix_column.isin(y.index)].drop_duplicates("matrix_column").set_index("matrix_column")
    multi = acc in ("GSE142401", "GSE118280", "GSE120862")
    sh["unit"] = sh.subject_id.astype(str) + ("|" + sh.group_arm.astype(str) if multi else "")
    if DS[acc][3] == "descriptive":
        sh["elapsed_h"] = np.where(sh.timepoint_label.astype(str).str.contains("before|baseline|Pre", case=False, regex=True), 0.0, 1.0)
    if acc == "GSE303983":
        sh.loc[sh.elapsed_h > 0, "elapsed_h"] = 2.25
    sh = sh[sh.elapsed_h <= 12]
    ok = sh.groupby("unit").elapsed_h.transform(lambda v: (v == 0).sum() == 1 and (v > 0).sum() >= 1 and v.nunique() == len(v))
    sh = sh[ok]
    return y.loc[sh.index], sh


def states(y, wf):
    g = [x for x in wf.genes if x in y.columns]
    w = wf.subset(g)
    L = P.projector(y[g], w, counts=False)
    b = y[g].to_numpy(float) @ L.T
    return b[:, 0] + 1j * b[:, 1]


class Data:
    def __init__(self, acc, placebo=False):
        self.acc = acc
        self.y, self.sh = load(acc)
        self.unit = pd.factorize(self.sh.unit)[0]
        self.el = self.sh.elapsed_h.to_numpy(float).copy()
        if placebo:
            rng = np.random.default_rng([SEED, 9, sum(map(ord, acc))])
            for u in np.unique(self.unit):
                i = np.flatnonzero(self.unit == u)
                self.el[i] = self.el[rng.permutation(i)]
        self.t0 = DS[acc][1]
        self.genes = [g for g in WF.genes if g in self.y.columns]
        self.idx = [np.flatnonzero(self.unit == u) for u in np.unique(self.unit)]
        self.base = np.array([i[self.el[i] == 0][0] for i in self.idx])[self.unit]
        self.later = self.el > 0
        self.c = np.exp(1j * W * (self.t0 + self.el)) - np.exp(1j * W * self.t0) if np.isfinite(self.t0) else np.zeros(len(self.el), complex)

    def sums(self, z, el=None):
        """(sum Re(delta conj c), sum |c|^2); `el`: permuted elapsed labels (rows keep their samples)."""
        if el is None:
            d = z - z[self.base]
            m = self.later
            return float(np.sum((d[m] * np.conj(self.c[m])).real)), float(np.sum(abs(self.c[m]) ** 2))
        c = np.exp(1j * W * (self.t0 + el)) - np.exp(1j * W * self.t0)
        zb = np.empty(len(z), complex)
        for i in self.idx:
            zb[i] = z[i][el[i] == 0][0]
        m = el > 0
        return float(np.sum(((z - zb)[m] * np.conj(c[m])).real)), float(np.sum(abs(c[m]) ** 2))

    def perm_el(self, rng):
        el = self.el.copy()
        for i in self.idx:
            el[i] = self.el[rng.permutation(i)]
        return el

    def pseudo(self, kind, n, tag, wf):
        lm = self.y.mean(0).astype(float)
        lm = lm[~lm.index.duplicated()]
        g = [x for x in wf.genes if x in self.y.columns]
        return P.pseudo_waveforms(lm, n, np.random.default_rng([SEED, tag, sum(map(ord, self.acc))]), wf.subset(g), kind)


def rho(parts):
    return sum(p[0] for p in parts) / sum(p[1] for p in parts)


def p_up(nul, obs):
    return float((1 + np.sum(np.asarray(nul) >= obs)) / (1 + len(nul)))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    out = np.empty(len(p))
    out[o] = np.minimum(np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]), 1)
    return out


def run(args):
    n_perm, n_sets = (300, 5) if args.smoke else (N_PERM, N_SETS)
    D = {a: Data(a, args.placebo) for a in DS}
    for a, d in D.items():
        say(a, DS[a][0], "units", len(d.idx), "samples", len(d.el), "clock genes", len(d.genes), "later points", sorted(set(d.el[d.later])))
    tested = [a for a in DS if DS[a][3] != "descriptive"]
    wf2 = WF.drop(["PER1", "NR1D1"])
    Z = {a: states(D[a].y, WF) for a in DS}
    Z2 = {a: states(D[a].y, wf2) for a in DS}
    rng = np.random.default_rng([SEED, 1])
    nul = {a: [] for a in tested}
    nul2 = {a: [] for a in tested}
    for _ in range(n_perm):
        for a in tested:
            el = D[a].perm_el(rng)
            nul[a].append(D[a].sums(Z[a], el))
            nul2[a].append(D[a].sums(Z2[a], el))
    fam_null = lambda ds, N: np.array([rho([N[a][i] for a in ds]) for i in range(n_perm)])  # noqa: E731
    power = []
    for k, ds in FAMILIES.items():
        v = fam_null(ds, nul)
        power.append(dict(family=k, datasets=len(ds), units=sum(len(D[a].idx) for a in ds), deltas=int(sum(D[a].later.sum() for a in ds)),
                          clock_genes_min=min(len(D[a].genes) for a in ds), null_mean=float(v.mean()), null_sd=float(v.std()), mde=MDE_K * float(v.std())))
    Pw = pd.DataFrame(power)
    Pw["scored"] = Pw.mde <= MAX_MDE
    Pw.to_csv(OUT / "power.csv", index=False)
    pd.DataFrame([dict(dataset=a, tissue=DS[a][0], role=DS[a][3], t0=DS[a][1], units=len(D[a].idx), samples=len(D[a].el), clock_genes=len(D[a].genes),
                       later_points=";".join(f"{x:.2f}" for x in sorted(set(D[a].el[D[a].later])))) for a in DS]).to_csv(OUT / "design.csv", index=False)
    pd.set_option("display.width", 250)
    print(Pw.round(4).to_string(), flush=True)
    if args.stage == "power":
        return

    comp, comp2, pz = {}, {}, {}
    for kind, tag in (("matched", 3), ("responsive", 4)):
        comp[kind] = {a: [D[a].sums(states(D[a].y, w)) for w in D[a].pseudo(kind, n_sets, tag, WF)] for a in tested}
        comp2[kind] = {a: [D[a].sums(states(D[a].y, w)) for w in D[a].pseudo(kind, n_sets, tag + 10, wf2)] for a in tested}
        say("companions", kind)
    tests = []
    for label, ZZ, NN, CC in (("all genes", Z, nul, comp), ("without PER1 and NR1D1", Z2, nul2, comp2)):
        for k, ds in FAMILIES.items():
            if label != "all genes" and k.startswith("H252.3"):
                continue
            obs = rho([D[a].sums(ZZ[a]) for a in ds])
            v = fam_null(ds, NN)
            row = dict(family=k if label == "all genes" else k.replace("H252.1", "H252.4a").replace("H252.2", "H252.4b"), genes=label, rho=obs,
                       null_mean=float(v.mean()), null_sd=float(v.std()), p=p_up(v, obs), mde=MDE_K * float(v.std()))
            for kind in CC:
                s = np.array([rho([CC[kind][a][j] for a in ds]) for j in range(n_sets)])
                row.update({f"{kind}_median": float(np.median(s)), f"{kind}_q95": float(np.quantile(s, .95)), f"p_{kind}": p_up(s, obs)})
            row["clock_specific"] = bool(obs > row["matched_q95"] and obs > row["responsive_q95"])
            tests.append(row)
    T = pd.DataFrame(tests)
    for fam in (["H252.3 muscle", "H252.3 adipose"], [x for x in T.family if x.startswith("H252.4")]):
        m = T.family.isin(fam)
        T.loc[m, "p_holm"] = holm(T.loc[m, "p"])

    rows, desc = [], []
    for a in tested:
        d = D[a]
        delta = Z[a] - Z[a][d.base]
        for e in sorted(set(d.el[d.later])):
            m = d.el == e
            mean = delta[m].mean()
            c = d.c[m][0]
            implied = ((np.angle(mean) - np.pi / 2) / W) % 24
            rows.append(dict(dataset=a, tissue=DS[a][0], elapsed_h=e, units=int(m.sum()), rho=float(np.sum((delta[m] * np.conj(c)).real) / (m.sum() * abs(c) ** 2)),
                             modulus_mean_delta=float(abs(mean)), predicted_modulus_at_rho1=float(abs(c)), implied_midpoint_h=float(implied),
                             protocol_midpoint_h=float((d.t0 + e / 2) % 24), rho_without_PER1_NR1D1=float(np.sum(((Z2[a] - Z2[a][d.base])[m] * np.conj(c)).real) / (m.sum() * abs(c) ** 2))))
        for shift in (-2, 2):
            c = np.exp(1j * W * (d.t0 + shift + d.el)) - np.exp(1j * W * (d.t0 + shift))
            desc.append(dict(item=f"rho with t0 {shift:+d} h", dataset=a, value=float(np.sum((delta[d.later] * np.conj(c[d.later])).real) / np.sum(abs(c[d.later]) ** 2))))
    for shift in (-2, 2):
        ds = FAMILIES["H252.2 seven datasets, assumed times"]
        num = den = 0
        for a in ds:
            d = D[a]
            c = np.exp(1j * W * (d.t0 + shift + d.el)) - np.exp(1j * W * (d.t0 + shift))
            dl = (Z[a] - Z[a][d.base])[d.later]
            num, den = num + np.sum((dl * np.conj(c[d.later])).real), den + np.sum(abs(c[d.later]) ** 2)
        desc.append(dict(item=f"H252.2 rho with t0 {shift:+d} h", dataset="pooled", value=float(num / den)))
    # rest vs exercise
    rng = np.random.default_rng([SEED, 6])
    for a, key in (("GSE252357", "CTRL"), ("GSE120862", "NE")):
        d = D[a]
        arm = d.sh.group_arm.astype(str).str.contains(key).to_numpy()
        delta = Z[a] - Z[a][d.base]
        times = [e for e in sorted(set(d.el[d.later])) if (arm & (d.el == e)).any() and (~arm & (d.el == e)).any()]
        stat = lambda r: float(np.mean([abs(delta[~r & (d.el == e)].mean() - delta[r & (d.el == e)].mean()) for e in times]))  # noqa: E731
        obs = stat(arm)
        ua = np.array([arm[i][0] for i in d.idx])
        nl = []
        for _ in range(2000 if not args.smoke else 100):
            if a == "GSE252357":
                pa = rng.permutation(ua)[d.unit]
            else:
                key2 = (d.sh.subject_id.astype(str) + d.sh.group_arm.astype(str).str.split("/").str[0]).to_numpy()
                flip = {k: rng.random() < 0.5 for k in np.unique(key2)}
                pa = arm ^ np.array([flip[k] for k in key2])
            nl.append(stat(pa))
        desc.append(dict(item="rest vs exercise: mean |difference of mean delta|", dataset=a, value=obs, null_mean=float(np.mean(nl)), p=p_up(nl, obs),
                         units_rest=int(ua.sum()), units_exercise=int((~ua).sum())))
        for lab, r in (("rest", arm), ("exercise", ~arm)):
            m = r & d.later
            desc.append(dict(item=f"rho, {lab} arm", dataset=a, value=float(np.sum((delta[m] * np.conj(d.c[m])).real) / np.sum(abs(d.c[m]) ** 2))))
    for a in [x for x in DS if DS[x][3] == "descriptive"]:
        d = D[a]
        mod = lambda z: float(abs((z - z[d.base])[d.later].mean()))  # noqa: E731
        obs = mod(Z[a])
        mean = (Z[a] - Z[a][d.base])[d.later].mean()
        row = dict(item="modulus of the mean delta", dataset=a, value=obs, units=len(d.idx), implied_midpoint_h=float(((np.angle(mean) - np.pi / 2) / W) % 24),
                   value_without_PER1_NR1D1=mod(Z2[a]))
        for kind, tag in (("matched", 3), ("responsive", 4)):
            s = np.array([mod(states(d.y, w)) for w in d.pseudo(kind, n_sets, tag, WF)])
            row.update({f"{kind}_median": float(np.median(s)), f"{kind}_q95": float(np.quantile(s, .95))})
        desc.append(row)
    T.to_csv(OUT / "score_tests.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "score_datasets.csv", index=False)
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    print(T.round(4).T.to_string(), pd.DataFrame(rows).round(3).to_string(), pd.DataFrame(desc).round(3).to_string(), sep="\n", flush=True)


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
