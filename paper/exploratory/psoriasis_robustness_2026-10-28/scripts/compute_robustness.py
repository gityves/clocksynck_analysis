"""Post hoc, descriptive: gene-set sensitivity of the psoriasis core / cue offsets (A54 held-out, A78, A81b H81.3).

Imports the existing estimator (figure_5/scripts/psoriasis/_shared, _cues; clocksync.paired) and changes nothing in it.
Pooled estimate = weighted resultant of per-series complex correlations (lesional vs uninvolved), weights 1/se^2 of
the A54 variance ratio, fixed at the headline weights. No permutations. Seed 2026.

    python compute_robustness.py
"""
from __future__ import annotations

import os
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
PSO = ROOT / "paper/results/figure_5/scripts/psoriasis"
OUT = Path(__file__).resolve().parents[1] / "output"
OLD = ROOT / "paper/results/figure_5/output/psoriasis"
sys.path.insert(0, str(PSO))
import numpy as np
import pandas as pd
import _shared as H
import _cues as Q
import clocksync as cs
from clocksync import stats as S
from clocksync.genes import canonical, ALIASES

PD, W, HOURS = H.PD, H.W, H.HOURS
SEED = 2026
A124 = ["FAM46B", "MTERFD3", "HIST2H2BE", "C12orf23", "SEPT4", "C10orf116", "CENPC1", "FAM134B"]


def wrap(x):
    return float(S.wrap(x))


class Series:
    """Paired tables of the gated series of one weight set."""
    def __init__(self, samples, w):
        self.names, self.pairs, self.cnt, self.w = [], {}, {}, []
        for k, g in samples.groupby("cohort", sort=False):
            if k in w:
                self.names.append(k)
                self.pairs[k] = H.split_pairs(g)
                self.cnt[k] = H.is_counts(g)
                self.w.append(w[k])
        self.w = np.array(self.w)

    def offset(self, wf_of):
        """Pooled offset (h); wf_of(cohort, case_columns) -> Waveform or callable."""
        cc = []
        for k in self.names:
            case, ref = self.pairs[k]
            wf = wf_of(k, case)
            cc.append(PD.ccorr(*PD.cohort_states(case, ref, self.cnt[k], wf)))
        return PD.pooled_offset(cc, self.w)


def sub_core(genes):
    return lambda k, case: H.WF.subset([g for g in genes if g in case.columns])


def sub_cue(wf, genes):
    return lambda k, case: H.restrict(wf, [g for g in genes if g in wf.genes and g in case.columns])


def summ(x):
    x = np.asarray(x, float)
    return dict(q025=np.quantile(x, .025), q50=np.quantile(x, .5), q975=np.quantile(x, .975))


def chunked_means(cohorts, chunk=4000):
    """Per-series mean log expression of the uninvolved samples of the paired patients, as the builders' `gene_means`
    (same rules, incl. the reference means of the no-level series) but reading the all-gene parquet in column chunks (the
    wide 88k-gene tables exhaust memory when read whole)."""
    import pyarrow.parquet as pq
    T = H.T
    bd, bh = H.builder("skin_paired_discovery"), H.builder("skin_paired_heldout")
    ref = pd.read_csv(bh.SRC / "GSE30999_nonlesional_means.csv", index_col=0).iloc[:, 0]
    out = {}
    for k in cohorts:
        src = bh.SRC if (bh.SRC / f"{k.split('|')[0]}.parquet").exists() else bd.SRC
        f = src / f"{k.split('|')[0]}.parquet"
        names = pq.ParquetFile(f).schema.names
        meta = [c for c in names if c in T.META]
        genes = [c for c in names if c not in T.META]
        m = pd.read_parquet(f, columns=meta)
        keep = ((m.disease == "psoriasis") & (m.status == "nonlesional")).to_numpy().copy()
        les = set(m.patient[(m.disease == "psoriasis") & (m.status == "lesional")])
        keep &= m.patient.isin(les).to_numpy()
        vt = m.value_type[keep].iloc[0]
        if src is bh.SRC and vt in bh.NO_LEVEL:
            gn = ["ARNTL" if g == "BMAL1" else g for g in genes] if "ARNTL" not in names else genes
            out[k] = ref[ref.index.isin(gn)]
            continue
        parts = []
        for i in range(0, len(genes), chunk):
            d = pd.read_parquet(f, columns=meta + genes[i:i + chunk])
            d = d.rename(columns={"BMAL1": "ARNTL"}) if "ARNTL" not in names else d
            non = d[keep].set_index("patient").sort_index()
            parts.append(T.mean_log_expression(non, T.META, vt == "counts"))
        out[k] = pd.concat(parts)
    return out


def main():
    import warnings
    wfs = Q.cue_waveforms(H.reference_sets()["systemic"])
    # gated cohorts and weights (discovery then held-out, as H.gated_samples)
    wgt = {}
    for which in ("discovery", "heldout"):
        p_ = pd.read_csv(OLD / f"a54_{which}_per_series.csv").set_index("cohort")
        wgt.update({k: 1 / p_.loc[k, "se"] ** 2 for k in p_.index if bool(p_.loc[k, "gated"])})
    means = chunked_means(wgt)
    cue_core = {"cue": wfs["cue"], "core": H.WF}
    # pseudo sets are drawn first (they need only the per-series mean expression) so that only the columns needed are
    # kept in memory (the all-gene tables of 15 series exceed 30 GB)
    draws, wanted = {}, set(H.GENES) | set(wfs["cue"].genes) | set(wfs["systemic"].genes)
    for kind in ("responsive", "matched"):
        for base_name, base in cue_core.items():
            if (kind, base_name) == ("responsive", "cue"):
                continue  # responsive pool (<= 27 present genes) is smaller than the 37-40 cue genes: sets would repeat genes, which the projector cannot take
            rg = np.random.default_rng([SEED, 1 if kind == "responsive" else 2, 1 if base_name == "cue" else 2])
            for k in wgt:
                lm = means[k]
                g = [x for x in base.genes if x in lm.index]
                wf = H.restrict(base, g) if base_name == "cue" else base.subset(g)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    sets = PD.pseudo_waveforms(lm, 200, rg, wf, kind)
                draws[(kind, base_name, k)] = sets
    meta_cols = None
    orig = H.T.read_all_genes
    allgenes = set().union(*[set(m.index) for m in means.values()])

    def slim(path, meta):
        """Read only the META and wanted columns of the parquet (BMAL1 renamed ARNTL as `read_all_genes`)."""
        import pyarrow.parquet as pq
        names = pq.ParquetFile(path).schema.names
        want = set(meta) | wanted | ({"BMAL1"} if "ARNTL" in wanted else set())
        d = pd.read_parquet(path, columns=[c for c in names if c in want])
        return d.rename(columns={"BMAL1": "ARNTL"}) if "ARNTL" not in names else d
    H.T.read_all_genes = slim
    samples, w = H.gated_samples(OLD)
    S15 = Series(samples, w)
    assert list(S15.names) == list(wgt), "cohort order"
    cue = wfs["cue"]
    core_genes = list(H.GENES)
    full_core = lambda k, case: H.WF.subset([g for g in core_genes if g in case.columns])
    full_cue = lambda k, case: H.restrict(cue, [g for g in cue.genes if g in case.columns])

    # ---- 1. reproduction -------------------------------------------------------------------------------
    core_o, cue_o = S15.offset(full_core), S15.offset(full_cue)
    sysm = S15.offset(lambda k, c: H.restrict(wfs["systemic"], [g for g in wfs["systemic"].genes if g in c.columns]))
    h3 = pd.read_csv(OLD / "a81b_h3_sets.csv").set_index("set").offset_h
    v78 = pd.read_csv(OLD / "a78_psoriasis_sets.csv").set_index("set").offset_h
    v3 = pd.read_csv(OLD / "a81b_h3_verdict.csv", header=None, index_col=0)[1]
    # held-out (A54), own weights
    ps = pd.read_csv(OLD / "a54_heldout_per_series.csv").set_index("cohort")
    ps = ps[ps.gated.astype(bool)]
    HO = Series(samples[samples.cohort.isin(ps.index)], {k: 1 / ps.loc[k, "se"] ** 2 for k in ps.index})
    ho_o = HO.offset(full_core)
    ho_old = float(pd.read_csv(OLD / "a54_heldout_pooled.csv", header=None, index_col=0)[1]["offset_h"])
    rows = [("A78/A81b core (15 series)", core_o, h3["core"]), ("A81b cue", cue_o, h3["cue"]),
            ("A78 systemic", sysm, v78["systemic"]),
            ("core minus cue", wrap(core_o - cue_o), float(v3["core_minus_cue_h"])),
            ("A54 held-out core (6 series)", ho_o, ho_old)]
    rep = pd.DataFrame(rows, columns=["estimate", "recomputed_h", "stored_h"])
    rep["abs_diff"] = (rep.recomputed_h - rep.stored_h).abs()
    rep.to_csv(OUT / "reproduction.csv", index=False)
    print(rep.to_string(index=False), flush=True)
    assert (rep.abs_diff < 1e-6).all(), "headline mismatch"
    head = dict(core=core_o, cue=cue_o, cc=wrap(core_o - cue_o), ho=ho_o)

    # ---- 2. leave-one-gene-out -------------------------------------------------------------------------
    rows = []
    for g in core_genes:
        gs = [x for x in core_genes if x != g]
        o, ho = S15.offset(sub_core(gs)), HO.offset(sub_core(gs))
        rows.append(dict(family="core", gene=g, core_h=o, cue_h=cue_o, core_minus_cue_h=wrap(o - cue_o), heldout_core_h=ho))
    present = sorted({g for k in S15.names for g in cue.genes if g in S15.pairs[k][0].columns})
    for g in present:
        gs = [x for x in present if x != g]
        o = S15.offset(sub_cue(cue, gs))
        rows.append(dict(family="cue", gene=g, core_h=core_o, cue_h=o, core_minus_cue_h=wrap(core_o - o), heldout_core_h=ho_o))
    lo = pd.DataFrame(rows)
    lo["d_core_h"] = [wrap(a - head["core"]) for a in lo.core_h]
    lo["d_cc_h"] = [wrap(a - head["cc"]) for a in lo.core_minus_cue_h]
    lo["d_heldout_h"] = [wrap(a - head["ho"]) for a in lo.heldout_core_h]
    lo["n_series_with_gene"] = [sum(g in S15.pairs[k][0].columns for k in S15.names) for g in lo.gene]
    lo.to_csv(OUT / "leave_one_gene_out.csv", index=False)
    rep_lo = []
    for fam, col, d in (("core", "core_h", "d_core_h"), ("core", "core_minus_cue_h", "d_cc_h"), ("core", "heldout_core_h", "d_heldout_h"),
                        ("cue", "cue_h", None), ("cue", "core_minus_cue_h", "d_cc_h")):
        x = lo[lo.family == fam]
        dd = x[d] if d else pd.Series([wrap(a - head["cue"]) for a in x.cue_h], index=x.index)
        infl = ";".join(f"{g}({v:+.2f})" for g, v in zip(x.gene[dd.abs() > .5], dd[dd.abs() > .5]))
        rep_lo.append(dict(family=fam, estimate=col, n_runs=len(x), headline=head[{"core_h": "core", "core_minus_cue_h": "cc",
                           "heldout_core_h": "ho", "cue_h": "cue"}[col]], min=x[col].min(), median=x[col].median(), max=x[col].max(),
                           max_abs_shift=dd.abs().max(), influential_gt_0p5h=infl))
    pd.DataFrame(rep_lo).to_csv(OUT / "leave_one_gene_out_summary.csv", index=False)
    print(pd.DataFrame(rep_lo).to_string(index=False), flush=True)

    # ---- 3. random 80% subsets -------------------------------------------------------------------------
    rng = np.random.default_rng(SEED)
    rows = []
    for j in range(500):
        k = int(round(.8 * len(present)))
        gs = list(rng.choice(present, k, replace=False))
        o = S15.offset(sub_cue(cue, gs))
        rows.append(dict(family="cue", draw=j, n_genes=k, core_h=core_o, cue_h=o, core_minus_cue_h=wrap(core_o - o)))
    for j in range(500):
        gs = list(rng.choice(core_genes, int(round(.8 * 14)), replace=False))
        o = S15.offset(sub_core(gs))
        rows.append(dict(family="core", draw=j, n_genes=len(gs), core_h=o, cue_h=cue_o, core_minus_cue_h=wrap(o - cue_o)))
    sb = pd.DataFrame(rows)
    sb.to_csv(OUT / "random_subsets.csv", index=False)
    rs = []
    for fam in ("cue", "core"):
        x = sb[sb.family == fam]
        for col, h in (("core_minus_cue_h", head["cc"]), ("core_h", head["core"])) + ((("cue_h", head["cue"]),) if fam == "cue" else ()):
            rs.append(dict(family=fam, estimate=col, headline=h, **summ(x[col]), share_same_sign=float((np.sign(x[col]) == np.sign(h)).mean()),
                           share_below_minus1h=float((x[col] < -1).mean()) if col != "cue_h" else np.nan))
    pd.DataFrame(rs).to_csv(OUT / "random_subsets_summary.csv", index=False)
    print(pd.DataFrame(rs).round(3).to_string(index=False), flush=True)

    # ---- 4. pseudo sets ----------------------------------------------------------------------------------
    # one series at a time (all-gene table of that file only), same table construction as H.discovery / H.heldout
    files = {}
    for bn in ("skin_paired_discovery", "skin_paired_heldout"):
        bb = H.builder(bn)
        files.update({f"{f.stem}|psoriasis": (bb.SRC, f.name) for f in sorted(Path(bb.SRC).glob("GSE*.parquet"))})
    cc = {(kd, bn): [[] for _ in range(200)] for kd in ("responsive", "matched") for bn in cue_core if (kd, bn) != ("responsive", "cue")}
    for k in S15.names:
        src, fn = files[k]
        need = set(H.GENES) | set(cue.genes)
        for (kd, bn, kk), sets in draws.items():
            if kk == k:
                for s_ in sets:
                    need |= set(s_.genes)
        wanted.clear()
        wanted.update(need)
        d = H.lesional_all_genes(src, [fn])
        d = d[d.cohort == k]
        case, ref = H.split_pairs(d)
        case, ref = (x.loc[:, ~x.columns.duplicated()] for x in (case, ref))  # duplicated symbols in a platform: first column
        cnt = H.is_counts(d)
        assert cnt == S15.cnt[k] and list(case.index) == list(S15.pairs[k][0].index) and list(ref.index) == list(S15.pairs[k][1].index)
        for (kd, bn), lst in cc.items():
            for j in range(200):
                wf = draws[(kd, bn, k)][j]
                assert all(g in case.columns for g in wf.genes)
                lst[j].append(PD.ccorr(*PD.cohort_states(case, ref, cnt, wf)))
        print("pseudo", k, flush=True)
        del d, case, ref
        import gc; gc.collect()
    rows, summ_rows = [], []
    for (kind, base_name), lst in cc.items():
        po = np.array([PD.pooled_offset(c, S15.w) for c in lst])
        d_core = np.array([wrap(core_o - p) for p in po])
        d_cue = np.array([wrap(cue_o - p) for p in po])
        for j in range(200):
            rows.append(dict(kind=kind, waveform=base_name, set=j, pseudo_offset_h=po[j], core_minus_pseudo_h=d_core[j]))
        summ_rows.append(dict(kind=kind, waveform=base_name, n_sets=200, **{"pseudo_" + a: b for a, b in summ(po).items()},
                              **{"core_minus_pseudo_" + a: b for a, b in summ(d_core).items()},
                              share_pseudo_le_core=float((po <= core_o).mean()), share_core_minus_pseudo_lt0=float((d_core < 0).mean()),
                              cue_minus_pseudo_q50=float(np.median(d_cue)),
                              mean_genes_per_set=float(np.mean([len(draws[(kind, base_name, k)][0].genes) for k in S15.names]))))
        print(summ_rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(OUT / "pseudo_sets.csv", index=False)
    pd.DataFrame(summ_rows).to_csv(OUT / "pseudo_sets_summary.csv", index=False)

    # ---- 5. platform names -------------------------------------------------------------------------------
    import pyarrow.parquet as pq
    cov = []
    cue_all = list(cue.genes)
    b1, b2 = H.builder("skin_paired_discovery"), H.builder("skin_paired_heldout")
    files = {f"{f.stem}|psoriasis": f for b in (b1, b2) for f in sorted(Path(b.SRC).glob("GSE*.parquet"))}
    for k in S15.names:
        cols = [c for c in pq.ParquetFile(files[k]).schema.names if c in allgenes or c == "BMAL1"]
        cset = set(cols)
        canon = {canonical(c) for c in cols}
        row = dict(series=k, n_pairs=len(S15.pairs[k][0]), platform_genes=len(cset),
                   core_present=sum(g in cset for g in core_genes),
                   core_present_canonical=sum(canonical(g) in canon for g in core_genes),
                   cue_present=sum(g in cset for g in cue_all), cue_present_canonical=sum(canonical(g) in canon for g in cue_all),
                   n_cue=len(cue_all),
                   a124_aliases_in_columns=";".join(a for a in A124 if a in cset),
                   a124_canonical_in_columns=";".join(a for a in A124 if canonical(a) in canon and canonical(a) != a))
        cov.append(row)
    pd.DataFrame(cov).to_csv(OUT / "platform_coverage.csv", index=False)
    al = pd.DataFrame([dict(alias=a, in_ALIASES=a in ALIASES, canonical=canonical(a), canonical_is_cue=canonical(a) in cue_all,
                            alias_is_cue=a in cue_all, canonical_is_core=canonical(a) in core_genes) for a in A124])
    al.to_csv(OUT / "a124_alias_check.csv", index=False)
    print(pd.DataFrame(cov).to_string(index=False), "\n", al.to_string(index=False))


if __name__ == "__main__":
    main()
