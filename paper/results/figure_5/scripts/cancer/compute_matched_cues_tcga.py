"""A84, TCGA arm: do tissue-matched GTEx cue sets reveal a core-versus-cue offset in tumours?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A84; H84, feasibility readout, secondary per cancer). Tumour minus adjacent
normal of the gated TCGA pairs (A52 units, weights 1 / se^2 of the cancer's variance ratio). Gene sets (A81b pools fixed: the
150 eligible genes of the 11-tissue selection table, split at the median lock into cue and locked candidates): for the tissue T
of each cancer (COAD colon = Transverse and Sigmoid pooled by complex mean, LIHC liver, BRCA breast, PRAD prostate, KIRC and KIRP
kidney, UCEC uterus, LUAD lung, THCA thyroid; BLCA has no GTEx tissue: the global A81b sets), candidates expressed in T with
|phase_T - global phase| < 3 h and p_fast,T < 0.2, the top 40 by amplitude in T, amplitude and phase from T; fewer than 15 genes:
the global A81b set restricted to genes expressed in T (`fallback`). Estimand: the pooled offset (h) of the core clock, the cue set
and the locked set (circular, weights), core minus cue (reference) and core minus locked, with identical per-patient tumour /
normal label swaps across sets (10,000; 1,000 per cancer), two-sided p. H84 (TCGA part): core minus matched cue, two-sided p
< 0.05 (the Holm over IBD and TCGA of the registration needs the IBD arm, not run here). Feasibility readout: weighted
cue |complex correlation|, global A81b set vs matched (matched sets are "more powerful" only if it rises).

    python compute_matched_cues_tcga.py [--subsample N] [--seed S]     (needs compute_tcga_paired.py run with the same --subsample;
        and figure_S07's compute_cue_reference.py for the selection tables)
Outputs (paper/results/figure_5/output/cancer/): a84_sets_info.csv, a84_set_<tissue>_<cue|locked>.csv, a84_tcga_summary.csv,
a84_tcga_per_cancer.csv, a84_tcga_feasibility.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a84_matched_cues.py (TCGA part; `a82_cues_diseases.tcga`).
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import ccorr, cohort_states, pooled_offset  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer"
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_out  # noqa: E402
S07 = group_out("gtex_cues")
W = 2 * np.pi / 24
N_SET, MIN_SET = 40, 15
TISSUE = {"colon": ["Colon - Transverse", "Colon - Sigmoid"], "liver": ["Liver"], "breast": ["Breast - Mammary Tissue"],
          "prostate": ["Prostate"], "kidney": ["Kidney - Cortex"], "uterus": ["Uterus"], "lung": ["Lung"],
          "thyroid": ["Thyroid"]}
CANCER = {"COAD": "colon", "LIHC": "liver", "BRCA": "breast", "PRAD": "prostate", "KIRC": "kidney", "KIRP": "kidney",
          "UCEC": "uterus", "LUAD": "lung", "THCA": "thyroid", "BLCA": None}
CLOCK = set(cs.gtex().genes) | {"BMAL1"}
CLOCK_OUTPUT = {"HLF", "NFIL3", "BHLHE40", "BHLHE41", "RORC", "NPAS2", "ARNTL2", "CRY2", "TSPAN4", "USP2", "LONRF3",
                "TEF", "NAMPT", "KLF10", "DDIT4L", "WEE1", "PPARA", "CIPC"}
EXCLUDE = CLOCK | CLOCK_OUTPUT  # local copy of figure_S07/_genes.py (same name `_shared` as this folder's module)


def per_gene(g: pd.DataFrame) -> pd.DataFrame:
    """Across tissues: number of tissues, Fisher combined p, weighted coherence Rw, lock (median over tissues with p_fast
    < 0.05), amplitude and acrophase (h) of the summed complex loading (local copy of figure_S07/_genes.py `per_gene`)."""
    from scipy import stats
    g = g.assign(c=g.amp_fast * np.exp(1j * W * g.acro_h))
    rows = {}
    for gene, d in g.groupby("gene"):
        rows[gene] = dict(nt=len(d), fisher=stats.chi2.sf(-2 * np.log(d.p_fast).sum(), 2 * len(d)),
                          Rw=abs(d.c.sum()) / d.c.abs().sum(), lock=d.loc[d.p_fast < 0.05, "lock"].median(),
                          amp=abs(d.c.sum()) / len(d), acro_h=np.angle(d.c.sum()) / W % 24)
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("gene")


def wf_of(x: pd.DataFrame, name: str) -> cs.Waveform:
    names = tuple(canonical(n) or n for n in x.index)
    return cs.Waveform(names, x.amp.to_numpy(float), x.acro_h.to_numpy(float) * W, None, None, name)


def pools() -> dict:
    """A81b candidate pools: eligible genes of the 11-tissue table split at the median lock."""
    agg = per_gene(pd.read_csv(S07 / "a81b_selection_genes.csv.gz"))
    el = agg[(agg.nt >= 8) & (agg.fisher < 1e-4) & (agg.Rw >= 0.6) & ~agg.index.isin(EXCLUDE)]
    med = el.lock.median()
    return {"cue": el[el.lock < med], "locked": el[el.lock >= med]}


def matched(g: pd.DataFrame, pool: dict, tissues, glob: pd.DataFrame):
    """Tissue-matched cue and locked sets (and a table of their sizes; colon = the two segments pooled by complex mean)."""
    t = g[g.tissue.isin(tissues)].assign(c=lambda d: d.amp_fast * np.exp(1j * W * d.acro_h))
    t = t.groupby("gene").agg(c=("c", "mean"), p=("p_fast", "min"))
    out, rows = {}, []
    for k, cand in pool.items():
        x = t.loc[t.index.intersection(cand.index)].copy()
        x["amp"], x["acro_h"] = x.c.abs(), np.angle(x.c) / W % 24
        dphi = np.abs(((x.acro_h - cand.loc[x.index, "acro_h"]) + 12) % 24 - 12)
        sel = x[(dphi < 3) & (x.p < 0.2)].nlargest(N_SET, "amp")
        fallback = len(sel) < MIN_SET
        if fallback:  # global A81b set restricted to genes expressed in T
            sel = glob[(glob.set == k) & glob.index.isin(x.index)][["amp", "acro_h"]]
        out[k] = sel[["amp", "acro_h"]]
        rows.append(dict(set=k, genes=len(sel), fallback=fallback))
    return out, pd.DataFrame(rows)


def _swap(b, f):
    return (np.where(f[:, None], b[2], b[0]), np.where(f, b[3], b[1]), np.where(f[:, None], b[0], b[2]),
            np.where(f, b[1], b[3]))


def compare(items, w, sets: list[dict], n_perm: int, seed: int, ref: str = "cue", other: str = "locked") -> tuple[dict, pd.DataFrame]:
    """Pooled offsets per gene set with identical label swaps; `sets[i]` = set name -> Waveform of item i
    (local: `paired.compare_sets` takes one waveform dict for all cohorts; here each cancer has its own sets)."""
    names = list(sets[0])
    states = {k: [] for k in names}
    used = {k: [] for k in names}
    for (cancer, ct, cn), wfs in zip(items, sets):
        for k in names:
            wf = wfs[k]
            g = [x for x in wf.genes if x in ct.columns and x in cn.columns]
            wf = wf.subset(g) if k == "core" else cs.Waveform(tuple(g), np.asarray(wf.amp)[[wf.genes.index(x) for x in g]],
                                                              np.asarray(wf.psi)[[wf.genes.index(x) for x in g]], None, None, wf.name)
            states[k].append(tuple(np.asarray(v) for v in cohort_states(ct, cn, True, wf)))
            used[k].append(len(g))
    pooled = lambda st: pooled_offset([ccorr(*b) for b in st], w)
    obs = {k: pooled(v) for k, v in states.items()}
    rng = np.random.default_rng(seed)
    null = {k: np.empty(n_perm) for k in names}
    for p in range(n_perm):
        flips = [rng.random(len(b[0])) < 0.5 for b in states[names[0]]]
        for k, st in states.items():
            null[k][p] = pooled([_swap(b, f) for b, f in zip(st, flips)])
    two = lambda n, d: float(min(1.0, 2 * min(cs.stats.perm_p(n, d, "less", tie=0.0), cs.stats.perm_p(-n, -d, "less", tie=0.0))))
    diff = float(cs.stats.wrap(obs["core"] - obs[ref]))
    nd = cs.stats.wrap(null["core"] - null[ref])
    d_o, nd_o = float(cs.stats.wrap(obs["core"] - obs[other])), cs.stats.wrap(null["core"] - null[other])
    rows = [dict(set=k, offset_h=obs[k], p_lower=cs.stats.perm_p(null[k], obs[k], "less", tie=0.0),
                 corr_abs=float(abs(np.sum(w * np.array([ccorr(*b) for b in states[k]])) / w.sum())),
                 genes_min=min(used[k]), genes_max=max(used[k])) for k in names]
    v = dict(core_minus_cue_h=diff, p_lower=cs.stats.perm_p(nd, diff, "less", tie=0.0), p_two=two(nd, diff),
             null_q025=float(np.quantile(nd, 0.025)), null_q975=float(np.quantile(nd, 0.975)),
             core_minus_locked_h=d_o, p_two_locked=two(nd_o, d_o), cancers=len(items), pairs=int(sum(len(i[1]) for i in items)))
    return v, pd.DataFrame(rows)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm, n_each = (100, 50) if run["quick"] else (10_000, 1000)
    g26 = pd.read_csv(S07 / "selection_genes_26.csv.gz")
    glob_ref = pd.read_csv(S07 / "a81b_cue_reference.csv", index_col=0)
    pool = pools()
    S, info = {}, []
    for k, tis in TISSUE.items():
        S[k], r = matched(g26, pool, tis, glob_ref)
        info.append(r.assign(tissue=k))
        for name in ("cue", "locked"):
            S[k][name].assign(tissue=k, set=name).to_csv(out / f"a84_set_{k}_{name}.csv")
    pd.concat(info).to_csv(out / "a84_sets_info.csv", index=False)
    core = cs.gtex()
    glob = {"core": core, **{k: wf_of(x, k) for k, x in glob_ref.groupby("set")}}
    sets_of = {k: {"core": core, **{n: wf_of(S[k][n], n) for n in ("cue", "locked")}} for k in S}

    pc = pd.read_csv(out / "a52_per_cancer.csv").set_index("cancer")
    samples = sh.load_tcga(args.subsample, args.seed)
    items, w, per_set = [], [], []
    for c in sh.CANCERS:
        if c not in pc.index or not bool(pc.loc[c, "gated"]):
            continue
        case, ref = sh.pair_frames(samples, c)
        ct, cn = sh.tcga_full_counts(c, case.index)
        ct, cn = ct.assign(depth=ct.sum(axis=1)), cn.assign(depth=cn.sum(axis=1))
        items.append((c, ct, cn))
        w.append(1 / pc.loc[c, "se"] ** 2)
        per_set.append(sets_of[CANCER[c]] if CANCER[c] else glob)
    w = np.asarray(w)
    v, sets = compare(items, w, per_set, n_perm, 84)
    v["H84_tcga"] = "offset" if v["p_two"] < 0.05 else "no offset"
    pd.Series(v).to_csv(out / "a84_tcga_summary.csv", header=False)
    sets.to_csv(out / "a84_tcga_sets.csv", index=False)
    _, gsets = compare(items, w, [glob] * len(items), 200 if run["quick"] else 2000, 82)  # global sets (A82 units)
    feas = pd.DataFrame([dict(readout="cue |corr|, global A81b set", value=float(gsets.set_index("set").loc["cue", "corr_abs"])),
                         dict(readout="cue |corr|, matched sets", value=float(sets.set_index("set").loc["cue", "corr_abs"]))])
    feas.to_csv(out / "a84_tcga_feasibility.csv", index=False)
    rows = []
    for i, it in enumerate(items):
        vi, si = compare([it], w[[i]], [per_set[i]], n_each, 84)
        rows.append(dict(cancer=it[0], tissue=CANCER[it[0]] or "global", **vi, **{f"{r.set}_offset_h": r.offset_h for r in si.itertuples()}))
    pd.DataFrame(rows).to_csv(out / "a84_tcga_per_cancer.csv", index=False)
    print(pd.read_csv(out / "a84_sets_info.csv").to_string(index=False), "\n", pd.Series(v).to_string(), "\n", sets.round(3).to_string(),
          "\n", feas.round(3).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
