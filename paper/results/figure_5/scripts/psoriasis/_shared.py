"""Helpers shared by psoriasis and psoriasis_details: the paired-skin tables, the A54 pipeline
(gate, variance ratio, pooled offset, gene checks, composition split) on top of `clocksync.paired`, the A56-A57
score and visit helpers, and the (temporary) readers of the legacy single-cell tables.

Method calls are clocksync's. Everything marked "local: not in clocksync" is a small helper for a step the package
has no public call for (listed in the folder READMEs).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS  # noqa: E402
sys.path.insert(0, str(DATASETS / "_bulk_common"))
import bulk_tables as T  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import stats as S  # noqa: E402
from common.paths import DATA, out_dir  # noqa: E402

PD = importlib.import_module("clocksync.paired")  # the module (`cs.paired` is the function of the same name)

HOURS = S.HOURS
W = 2 * np.pi / 24
WF = cs.gtex()
GENES = list(WF.genes)
LEVELS = ("nonlesional", "lesional")
FOLDER = "psoriasis"
FOLDER_S = "psoriasis_details"

DISC_ORDER = ["GSE30999", "GSE117468", "GSE201827", "GSE136757", "GSE85034", "GSE252029", "GSE197056", "GSE268958",
              "GSE121212"]  # legacy a54_skin.DISC: the order fixes the random stream
KER = ["KRT1", "KRT10", "KRT16", "KRT6A", "LOR", "FLG"]
DERM = ["COL1A1", "COL1A2", "DCN"]
IL17 = ["IL17A", "IL36G", "DEFB4A", "S100A7A", "S100A9", "SERPINB4", "KRT16", "PI3", "LCN2"]
IFNG = ["CXCL9", "CXCL10", "IFNG", "STAT1"]
PROLIF = "MKI67 TOP2A CENPF NUSAP1 CCNB1 CDK1 BIRC5 TYMS MCM2 PCNA".split()
BASAL = "KRT5 KRT14 KRT15 COL17A1 TP63 ITGA6 ITGB4 CXCL14".split()
DIFF = "KRT1 KRT10 KRT2 DSG1 DSC1 SBSN KRTDAP LORICRIN LOR FLG IVL CNFN LCE1B LCE3D SPRR1A SPRR1B SPRR2A".split()


def builder(name: str):
    """Module `paper/common/datasets/<name>/build.py`."""
    f = DATASETS / name / "build.py"
    spec = importlib.util.spec_from_file_location(f"build_{name}", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---- tables -----------------------------------------------------------------------------------------------------

def thin(samples: pd.DataFrame, n: int | None, seed: int = 0, cohort: str = "cohort", unit: str = "unit") -> pd.DataFrame:
    """--subsample: keep n random units per cohort (all samples of a unit)."""
    if n is None:
        return samples
    rng = np.random.default_rng(seed)
    keep = []
    for k, g in samples.groupby(cohort, sort=False):
        u = sorted(g[unit].unique())
        keep += [(k, x) for x in rng.permutation(u)[:n]]
    key = pd.MultiIndex.from_tuples(keep)
    return samples[pd.MultiIndex.from_frame(samples[[cohort, unit]]).isin(key)].reset_index(drop=True)


def order_cohorts(samples: pd.DataFrame, order: list[str], col: str = "cohort") -> pd.DataFrame:
    """Rows of the listed cohorts in the given order (legacy processing order: it fixes the random stream)."""
    rank = {k: i for i, k in enumerate(order)}
    t = samples[samples[col].isin(rank)]
    return t.iloc[np.argsort(t[col].map(rank).to_numpy(), kind="stable")].reset_index(drop=True)


def lesional_all_genes(folder: Path, files=None, by_disease: bool = True) -> pd.DataFrame:
    """local: not in the builders. As `bulk_tables.lesional_series` but with every gene column, which the
    expression-matched pseudo-clock control of `clocksync.paired` needs (the builders return only the 14 clock genes
    and `gene_means()`; the pseudo gene sets must be columns of the sample table)."""
    files = sorted(Path(folder).glob("GSE*.parquet")) if files is None else [Path(folder) / f for f in files]
    out = []
    for f in files:
        d = T.read_all_genes(f, T.META)
        d.insert(0, "series", f.stem.split("_")[0])
        d["cohort"] = d.series + "|" + d.disease.astype(str) if by_disease else d.series
        out.append(d.rename(columns={"patient": "unit", "status": "condition"}).assign(
            unit=lambda x: x["unit"].astype(str)))
    return pd.concat(out, ignore_index=True)


def discovery(subsample=None, seed=0, all_genes=True):
    """(samples, means) of the nine A54 discovery psoriasis series, in the legacy order. `all_genes`: every gene
    column (needed with `pseudo_means`); False: the builder's 14-gene table."""
    b = builder("skin_paired_discovery")
    if all_genes:
        samples = lesional_all_genes(b.SRC)
    else:
        samples, _ = b.build()
    t = b.psoriasis_discovery(samples)
    t = order_cohorts(t, [f"{g}|psoriasis" for g in DISC_ORDER])
    return thin(t, subsample, seed), b.gene_means()


def heldout(subsample=None, seed=0, all_genes=True):
    b = builder("skin_paired_heldout")
    if all_genes:
        samples = lesional_all_genes(b.SRC)
        samples = samples[samples.disease == "psoriasis"].reset_index(drop=True)
    else:
        samples, _ = b.build()
    return thin(samples, subsample, seed), b.gene_means()


def split_pairs(x: pd.DataFrame, levels=LEVELS, unit="unit", condition="condition"):
    """local: not in clocksync (the private `paired._split`): case and reference rows aligned on unit."""
    case = x[x[condition].astype(str) == levels[1]].drop_duplicates(unit).set_index(unit)
    ref = x[x[condition].astype(str) == levels[0]].drop_duplicates(unit).set_index(unit)
    u = sorted(set(case.index) & set(ref.index))
    # genes absent from a series are all-NaN columns of the union table: drop them (legacy: per-series files)
    return case.loc[u].dropna(axis=1, how="all"), ref.loc[u].dropna(axis=1, how="all")


def is_counts(x: pd.DataFrame) -> bool:
    return bool((x.value_type == "counts").all())


def logexpr(x: pd.DataFrame, genes, counts: bool) -> np.ndarray:
    return PD.log_expr(x[list(genes)].astype(float), x.depth.to_numpy(float) if counts else None).to_numpy()


def states(case, ref, counts, waveform=WF):
    return PD.cohort_states(case, ref, counts, waveform)


def series_label(cohort: str) -> str:
    return cohort.split("|")[0]


# ---- A54 pipeline -----------------------------------------------------------------------------------------------

def gene_shift(case, ref, counts, g) -> complex:
    """kappa_L / kappa_N of gene g on the uninvolved host built without g (A54 H54.3b)."""
    wf = WF.subset([x for x in GENES if x != g])
    b = states(case, ref, counts, wf)
    h = S.complex_state(b[2])
    out = []
    for x in (case, ref):
        y = logexpr(x, [g], counts)[:, 0]
        out.append(np.sum((y - y.mean()) * np.conj(h)) / np.sum(abs(h) ** 2))
    return out[0] / out[1] if out[0] != 0 and out[1] != 0 else np.nan


def comp_score(case, ref, counts, pos=KER, neg=DERM) -> np.ndarray:
    """Keratinocyte-minus-dermal score, lesional minus uninvolved, per pair (A54 composition check)."""
    genes = [g for g in pos + neg if g in case.columns]
    y = np.vstack([logexpr(case, genes, counts), logexpr(ref, genes, counts)])
    z = pd.DataFrame((y - y.mean(0)) / y.std(0), columns=genes)
    s = z[[g for g in pos if g in genes]].mean(1) - z[[g for g in neg if g in genes]].mean(1)
    n = len(case)
    return s.iloc[:n].to_numpy() - s.iloc[n:].to_numpy()


def per_series_table(res) -> pd.DataFrame:
    d = res.cohorts.copy()
    d.insert(0, "series", d.cohort.map(series_label))
    d = d.rename(columns={"coherence": "ccorr_abs"})
    return d


def analyse(samples, means, tag, out: Path, seed: int, hyp: str, pos=KER, neg=DERM, n_perm=10000, n_boot=1000,
            min_pairs=10, genes_checks=True):
    """A54 pipeline (also A55, A59): gate on the uninvolved samples, variance ratio r with pseudo-clock control,
    pooled over gated series (swap p), pooled offset (swap p), leave-one-gene-out offsets, per-gene shifts and the
    composition split. Writes <tag>_per_series.csv, <tag>_pooled.csv and <tag>_genes.csv to `out`."""
    res = cs.paired(samples, unit="unit", condition="condition", levels=LEVELS, cohort="cohort", counts="value_type",
                    pseudo_means=means, n_perm=n_perm, n_boot=n_boot, min_pairs=min_pairs, seed=seed)
    d = per_series_table(res)
    d.to_csv(out / f"{tag}_per_series.csv", index=False)
    p = res.pooled
    if not p.get("n_cohorts", 0):
        raise ValueError("no gated series")
    keys = p.cohorts.split(";")
    use = res.cohorts.set_index("cohort").loc[keys]
    w = 1 / use.se.to_numpy() ** 2
    v = dict(series=";".join(series_label(k) for k in keys), n_pairs=int(p.n_pairs), r=p.fixed, se_fixed=p.se_fixed,
             r_re=p.random, lo_re=p.lo, hi_re=p.hi, tau2=p.tau2, k=int(p.k), p_swap_upper=p.p_swap_upper,
             pseudo_q95=p.get("pseudo_q95", np.nan), pseudo_median=p.get("pseudo_median", np.nan),
             offset_h=p.offset_h, p_offset_lower=p.p_offset_lower, p_offset_two_sided=p.p_offset)
    v[f"{hyp}_1"] = "PASS" if (p.fixed > 0 and p.p_swap_upper < 0.05 and p.fixed > v["pseudo_q95"]) else "FAIL"
    v[f"{hyp}_2"] = "PASS" if (p.offset_h < 0 and p.p_offset_lower < 0.05) else "FAIL"
    pairs = {k: split_pairs(samples[samples.cohort == k]) for k in keys}
    cnt = {k: is_counts(samples[samples.cohort == k]) for k in keys}
    if genes_checks:
        logo, shifts = {}, {}
        for g in GENES:
            wf = WF.subset([x for x in GENES if x != g])
            logo[g] = PD.pooled_offset([PD.ccorr(*states(*pairs[k], cnt[k], wf)) for k in keys], w)
            z = np.asarray([gene_shift(*pairs[k], cnt[k], g) for k in keys])
            ok = np.isfinite(z)
            shifts[g] = float(np.angle(np.sum(w[ok] * z[ok] / np.maximum(abs(z[ok]), 1e-12))) * HOURS)
        gl = pd.DataFrame(dict(logo_offset_h=logo, gene_shift_h=shifts)).rename_axis("gene")
        gl.to_csv(out / f"{tag}_genes.csv")
        v["logo_max_h"] = float(gl.logo_offset_h.max())
        v["n_genes_negative"] = int((gl.gene_shift_h < 0).sum())
        v[f"{hyp}_3"] = "PASS" if (v["logo_max_h"] < 0 and v["n_genes_negative"] >= 10) else "FAIL"
    lo, hi = [], []
    for k in keys:
        b = res.data[k]
        sd = comp_score(*pairs[k], cnt[k], pos, neg)
        m = sd <= np.median(sd)
        for sel, acc in ((m, lo), (~m, hi)):
            acc.append(PD.ccorr(*(x[sel] for x in b)))
    v["offset_low_composition_diff_h"] = PD.pooled_offset(lo, w)
    v["offset_high_composition_diff_h"] = PD.pooled_offset(hi, w)
    v["amp_lower_bound_wmean_log"] = float(np.sum(w * np.log(use.x_amp_lo.to_numpy())) / w.sum())
    pd.Series(v).to_csv(out / f"{tag}_pooled.csv", header=False)
    return res, d, v


# ---- weights from the A54 series (A56-A57, A78-A82) --------------------------------------------------------------

def gated_weights(per_series: pd.DataFrame, key="series") -> pd.Series:
    """1 / se^2 of the gated series (the A54 weights); index = series."""
    g = per_series[per_series.gated.astype(bool)].set_index(key)
    return 1 / g.se ** 2


def psoriasis_items(out: Path, subsample=None, seed=0, quick=False):
    """Gated A54 series (discovery and held-out): ({series: (case, ref, counts)}, weights). Weights are read from
    the A54 outputs in `out` (written by compute_replication.py; documented hand-over file)."""
    items, w = {}, {}
    for which, load in (("discovery", discovery), ("heldout", heldout)):
        ps = pd.read_csv(out / f"a54_{which}_per_series.csv").set_index("series")
        samples, _ = load(subsample, seed, all_genes=True)
        for k, g in samples.groupby("cohort", sort=False):
            s = series_label(k)
            if s in ps.index and bool(ps.loc[s, "gated"]):
                case, ref = split_pairs(g)
                items[s] = (case, ref, is_counts(g))
                w[s] = 1 / ps.loc[s, "se"] ** 2
    return items, pd.Series(w)


# ---- treatment visits ---------------------------------------------------------------------------------------------

def treatment_samples(file: str, all_genes: bool = False) -> pd.DataFrame:
    """local: not in the builder. The `skin_paired_treatment` builder's table of one file with the `pasi` column
    (`extra=("pasi",)` of `bulk_tables.lesional_series`); `all_genes`: every gene column (the A57 / A73 / A74 scores
    use genes outside the clock panel)."""
    src = DATA / "skin_lesional" / "processed_treatment"
    if all_genes:
        t = lesional_all_genes(src, [file])
    else:
        t = T.lesional_series(src, [file], extra=("pasi",))
        t["cohort"] = Path(file).stem + "|" + t.disease.astype(str)
    return t


def visits(file: str, treated: bool | None = None, drug: str | None = None, all_genes: bool = False):
    """{(visit, status): rows aligned on patient} of the patients with all four samples, and the counts flag."""
    d = treatment_samples(file, all_genes)
    if treated is not None:
        pl = d.arm.astype(str).str.lower().str.contains("placebo")
        d = d[~pl] if treated else d[pl]
    if drug:
        d = d[d.treatment.astype(str).str.contains(drug, case=False)]
    t = {(v, s): d[(d.visit == v) & (d.condition == s)].drop_duplicates("unit").set_index("unit")
         for v in ("baseline", "post") for s in LEVELS[::-1]}
    pts = sorted(set.intersection(*[set(x.index) for x in t.values()]))
    return {k: x.loc[pts] for k, x in t.items()}, bool(d.value_type.iloc[0] == "counts")


def visit_states(t, counts):
    """(baseline, post) states of lesional vs uninvolved pairs: tuples (bA, vA, bB, vB)."""
    b = states(t[("baseline", "lesional")], t[("baseline", "nonlesional")], counts)
    p = states(t[("post", "lesional")], t[("post", "nonlesional")], counts)
    return b, p


def wrap1(h) -> float:
    return float(S.wrap(h))


def score_diff(case, ref, counts, genes) -> np.ndarray:
    """z-scored mean of `genes` over the 2n samples, lesional minus uninvolved per pair (A56)."""
    g = [x for x in genes if x in case.columns]
    y = np.vstack([logexpr(case, g, counts), logexpr(ref, g, counts)])
    sd = y.std(0)
    z = (y - y.mean(0)) / np.where(sd > 0, sd, 1)
    s = z.mean(1)
    n = len(case)
    return s[:n] - s[n:]


def sub(b, m):
    return tuple(x[m] for x in b)


def offset_h(b) -> float:
    return float(np.angle(PD.ccorr(*b)) * HOURS)


def parser(doc: str):
    """`common.cli.parser` plus `--n-perm` / `--n-boot` (override the pre-registered resampling sizes, for checks)."""
    from common import cli
    ap = cli.parser(doc)
    ap.add_argument("--n-perm", type=int, default=None, help="permutations (default: as pre-registered)")
    ap.add_argument("--n-boot", type=int, default=None, help="bootstrap draws (default: as pre-registered)")
    return ap


def sizes(args, run: dict, n_perm: int = 10000, n_boot: int = 1000) -> dict:
    """Resampling sizes: pre-registered, 50 / 20 with --subsample, or the explicit overrides."""
    q = run["quick"]
    return dict(n_perm=args.n_perm or (50 if q else n_perm), n_boot=args.n_boot or (20 if q else n_boot))


# ---- timed healthy skin (GSE205155): per-gene cosinor and the A78 reference sets ----------------------------------

def timed_skin():
    """(samples, genes x samples matrix of log2 intensities collapsed on the canonical symbol) of GSE205155.
    Uses `_read` of the builder (the builder's public `build()` keeps only the 14 clock genes as pseudo-counts)."""
    from clocksync.genes import canonical
    s, m, sym = builder("skin_timed_gse205155")._read()
    c = sym.fillna("").map(canonical)
    keep = c.to_numpy() != ""
    return s, m[keep].groupby(c[keep].to_numpy()).mean()


def cosinor_f(y: np.ndarray, t: np.ndarray, subj: np.ndarray):
    """local: not in clocksync (`clocksync.timed.cosinor` has no test). y: samples x genes. Per gene, cosinor on clock time
    with subject intercepts: (acrophase h, amplitude, F-test p of the 24 h term)."""
    from scipy import stats as st
    Sd = pd.get_dummies(subj).to_numpy(float)
    X0, X1 = Sd, np.column_stack([Sd, np.cos(W * t), np.sin(W * t)])
    b1, *_ = np.linalg.lstsq(X1, y, rcond=None)
    r0 = y - X0 @ np.linalg.lstsq(X0, y, rcond=None)[0]
    r1 = y - X1 @ b1
    ss0, ss1 = (r0 ** 2).sum(0), (r1 ** 2).sum(0)
    df1 = len(y) - X1.shape[1]
    F = ((ss0 - ss1) / 2) / np.maximum(ss1 / df1, 1e-12)
    a, b = b1[-2], b1[-1]
    return np.mod(np.arctan2(b, a) / W, 24), np.hypot(a, b), st.f.sf(F, 2, df1)


def per_layer(s: pd.DataFrame, m: pd.DataFrame, genes) -> dict:
    out = {}
    for layer in ("epidermis", "dermis"):
        ss = s[s.layer == layer]
        out[layer] = cosinor_f(m.loc[genes, ss.gsm].to_numpy().T, ss.clock_h.to_numpy(float), ss.subject.to_numpy())
    return out


def reference_sets(out: Path | None = None, n_set: int = 40) -> dict:
    """A78 reference planes from the timed healthy epidermis: systemic (rhythmic in both layers p < 0.01, in phase
    within 0.5 h, top `n_set` epidermal amplitude) and epidermal-output (rhythmic in the epidermis p < 0.001 only,
    p_dermis > 0.05, top `n_set`) gene sets as waveforms (amplitude x ln 2, epidermal acrophase), plus the core clock."""
    s, m = timed_skin()
    genes = [g for g in m.index if g not in GENES]
    r = per_layer(s, m, genes)
    d = pd.DataFrame({"acro_e": r["epidermis"][0], "amp_e": r["epidermis"][1], "p_e": r["epidermis"][2],
                      "acro_d": r["dermis"][0], "p_d": r["dermis"][2]}, index=genes)
    d["shift"] = np.asarray(S.wrap(d.acro_d - d.acro_e))
    sysm = d[(d.p_e < 0.01) & (d.p_d < 0.01) & (d["shift"].abs() < 0.5)].nlargest(n_set, "amp_e")
    outp = d[(d.p_e < 0.001) & (d.p_d > 0.05) & ~d.index.isin(sysm.index)].nlargest(n_set, "amp_e")
    wf = lambda x, name: cs.Waveform(tuple(x.index), x.amp_e.to_numpy() * np.log(2), x.acro_e.to_numpy() * W, None,
                                     None, name)
    if out is not None:
        pd.concat([sysm.assign(set="systemic"), outp.assign(set="epidermal_output")]).to_csv(out / "a78_gene_sets.csv")
    return {"core": WF, "systemic": wf(sysm, "systemic"), "epidermal_output": wf(outp, "epidermal_output")}


def gated_samples(out: Path, subsample=None, seed=0):
    """Gated A54 psoriasis series (discovery then held-out) as one sample table with every gene column, and the
    cohort -> 1 / se^2 weights (A54 outputs in `out`, written by compute_replication.py)."""
    parts, w = [], {}
    for which, load in (("discovery", discovery), ("heldout", heldout)):
        ps = pd.read_csv(out / f"a54_{which}_per_series.csv").set_index("cohort")
        samples, _ = load(subsample, seed, all_genes=True)
        for k, g in samples.groupby("cohort", sort=False):
            if k in ps.index and bool(ps.loc[k, "gated"]):
                parts.append(g)
                w[k] = 1 / ps.loc[k, "se"] ** 2
    return pd.concat(parts, ignore_index=True), w


def restrict(wf, genes):
    idx = [wf.genes.index(g) for g in genes]
    return cs.Waveform(tuple(genes), np.asarray(wf.amp)[idx], np.asarray(wf.psi)[idx], None, None, wf.name)


def compare_cue_sets(samples: pd.DataFrame, weights: dict, waveforms: dict, reference: str, n_perm: int, seed: int,
                     cohort: str = "cohort", counts: str = "value_type"):
    """local: not in clocksync (`paired.compare_sets` with absent genes dropped per cohort; `waveforms` may be a callable
    cohort -> {set: Waveform} for tissue-matched sets). `compare_sets` keeps a gene
    whose column exists in the table; in a union table of series a gene absent from one series is an all-NaN column
    there, which makes the cohort's state NaN. Same estimator and the same swap stream as `compare_sets`: pooled
    offsets per gene set with identical per-pair label swaps across sets, set minus `reference` (wrapped)."""
    groups = [(n, x) for n, x in samples.groupby(cohort, sort=False) if n in weights]
    levels = LEVELS if "lesional" in set(samples.condition) else ("normal", "tumour")  # skin / gut or TCGA pairs
    sets_of = waveforms if callable(waveforms) else (lambda c: waveforms)  # a callable gives a set dict per cohort (A84)
    names = list(sets_of(groups[0][0]))
    w = np.array([weights[n] for n, _ in groups], float)
    states, used = {k: [] for k in names}, {k: [] for k in names}
    for name, x in groups:
        case, ref = split_pairs(x, levels)
        for k in names:
            wk = sets_of(name)[k]
            g = [y for y in wk.genes if y in case.columns and y in ref.columns]
            wf = wk.subset(g) if k == "core" else restrict(wk, g)
            states[k].append(tuple(np.asarray(v) for v in PD.cohort_states(case, ref, is_counts(x), wf)))
            used[k].append(len(g))
    pooled = lambda st: PD.pooled_offset([PD.ccorr(*b) for b in st], w)
    swp = lambda b, f: (np.where(f[:, None], b[2], b[0]), np.where(f, b[3], b[1]), np.where(f[:, None], b[0], b[2]),
                        np.where(f, b[1], b[3]))
    obs = {k: pooled(v) for k, v in states.items()}
    rng = np.random.default_rng(seed)
    null = {k: np.empty(n_perm) for k in names}
    for p in range(n_perm):
        flips = [rng.random(len(b[0])) < 0.5 for b in states[names[0]]]
        for k, st in states.items():
            null[k][p] = pooled([swp(b, f) for b, f in zip(st, flips)])
    sets = pd.DataFrame([dict(set=k, offset_h=obs[k], p_lower=S.perm_p(null[k], obs[k], "less", tie=0.0),
                              coherence=float(abs(np.sum(w * np.array([PD.ccorr(*b) for b in states[k]])) / w.sum())),
                              genes_min=min(used[k]), genes_max=max(used[k])) for k in names])
    rows = []
    for k in names:
        if k == reference:
            continue
        d, nd = float(S.wrap(obs[k] - obs[reference])), S.wrap(null[k] - null[reference])
        pl = S.perm_p(nd, d, "less", tie=0.0)
        rows.append(dict(set=k, reference=reference, difference_h=d, p_lower=pl,
                         p_two=float(min(1.0, 2 * min(pl, S.perm_p(-nd, -d, "less", tie=0.0)))),
                         null_q025=float(np.quantile(nd, 0.025)), null_q975=float(np.quantile(nd, 0.975)),
                         cohorts=len(groups), pairs=int(sum(len(b[0]) for b in states[names[0]]))))
    return sets, pd.DataFrame(rows)
