"""A254: is a clock delay of the lesion against the same person's unaffected tissue general outside the skin?

Registered in paper/docs/registered_experiments.md (A254). Estimator as A183 (`clocksync.paired` building blocks): both
samples of a dataset on one projector, offset = arg of the complex correlation of the paired states (h; negative = lesion
delayed), r = log ratio of the clock-plane variance (lesion / non-lesion). Null: labels swapped per person.

Stages
  power  pairs, gate, coherence flag, swap-null SD and MDE of each offset and of the pooled one -> output/power.csv
  score  (after the power note is dated) everything registered -> output/score_{datasets,pooled,descriptive}.csv
  `--placebo`: labels swapped at random per person once, before anything else (code test; output/placebo/).

    python compute_a254.py --stage power|score [--smoke] [--placebo]

Choices not fixed by the registration
- Input: `output/prep/<GSE>_{expr.parquet,sheet.csv}` (`prep_lesions_tissues.py`). Several samples of a person and condition
  are combined: counts summed, array log values averaged.
- A clock gene absent from a dataset is dropped from the waveform there.
- Companion sets are drawn from the non-lesional samples' mean log expression (counts: log mean fraction, genes with a mean
  count >= 1); set j is pooled across datasets with the same weights as the clock.
- Two-sided swap p: share of |null| >= |observed|; one-sided (delay): share of null <= observed. Offsets live in (-12, 12] h.
- Pooled r: inverse-variance (bootstrap SE, 1,000) fixed effect; its null pools the per-dataset swap nulls with the same weights.
- GSE39491: the lesion against each control tissue separately (descriptive).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
PREP = HERE.parent / "output" / "prep"
SEED, N_PERM, N_SETS, N_BOOT, MDE_K = 254, 10_000, 200, 1000, 2.49
WF = cs.gtex()
HOURS = 24 / (2 * np.pi)
TESTED = {"GSE16134": "periodontitis", "GSE57218": "osteoarthritis", "GSE43292": "carotid plaque", "GSE192352": "uterine fibroid",
          "GSE151243": "hidradenitis"}
DESC = {"GSE10334": ("periodontitis (second series)", "nonlesion"), "GSE39491|squamous": ("Barrett vs squamous oesophagus", "nonlesion_squamous"),
        "GSE39491|cardia": ("Barrett vs gastric cardia", "nonlesion_cardia")}


def say(*a):
    import time
    print(time.strftime("%F %T"), *a, flush=True)


class Pair:
    def __init__(self, key, ref_label="nonlesion", placebo=False):
        acc = key.split("|")[0]
        x = pd.read_parquet(PREP / f"{acc}_expr.parquet")
        s = pd.read_csv(PREP / f"{acc}_sheet.csv", dtype={"person": str}).set_index("sample")
        s = s[s.condition.isin(["lesion", ref_label])]
        x = x.loc[s.index]
        self.counts = bool((s.kind == "counts").all())
        agg = x.groupby([s.person.to_numpy(), s.condition.to_numpy()]).sum() if self.counts else x.groupby([s.person.to_numpy(), s.condition.to_numpy()]).mean()
        case, ref = agg.xs("lesion", level=1), agg.xs(ref_label, level=1)
        u = sorted(set(case.index) & set(ref.index))
        case, ref = case.loc[u].astype(float), ref.loc[u].astype(float)
        if placebo:
            sw = np.random.default_rng([SEED, 9, sum(map(ord, key))]).random(len(u)) < 0.5
            a, b = case.to_numpy().copy(), ref.to_numpy().copy()
            a[sw], b[sw] = ref.to_numpy()[sw], case.to_numpy()[sw]
            case, ref = pd.DataFrame(a, index=case.index, columns=case.columns), pd.DataFrame(b, index=ref.index, columns=ref.columns)
        if self.counts:
            case["depth"], ref["depth"] = case.sum(1), ref.sum(1)
        self.key, self.case, self.ref, self.n = key, case, ref, len(u)
        self.wf = WF.subset([g for g in WF.genes if g in x.columns])
        yr = P.log_expr(ref[list(self.wf.genes)], ref.depth.to_numpy(float) if self.counts else None).to_numpy()
        self.gate = P.reference_gate(yr, self.wf, 200, np.random.default_rng([SEED, 1, sum(map(ord, key))]))
        genes = [c for c in ref.columns if c != "depth"]
        if self.counts:
            self.lm = np.log(ref[genes].div(ref.depth, axis=0).mean() + 1e-12)[ref[genes].mean() >= 1]
        else:
            self.lm = (ref[genes].mean() * np.log(2.0)).dropna()
        self.lm = self.lm[~self.lm.index.duplicated()]

    def states(self, wf=None, drop=None):
        return P.cohort_states(self.case, self.ref, self.counts, wf or self.wf, drop_genes=drop)

    def pseudo(self, kind, n, tag):
        return P.pseudo_waveforms(self.lm, n, np.random.default_rng([SEED, tag, sum(map(ord, self.key))]), self.wf, kind)


def p_two(nul, obs):
    return float((1 + np.sum(np.abs(nul) >= abs(obs))) / (1 + len(nul)))


def p_low(nul, obs):
    return float((1 + np.sum(np.asarray(nul) <= obs)) / (1 + len(nul)))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    out = np.empty(len(p))
    out[o] = np.minimum(np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]), 1)
    return out


def run(args):
    n_perm, n_sets, n_boot = (300, 5, 100) if args.smoke else (N_PERM, N_SETS, N_BOOT)
    D = {k: Pair(k, placebo=args.placebo) for k in TESTED}
    rows = []
    B, nul_o, nul_r = {}, {}, {}
    for k, d in D.items():
        b = B[k] = d.states()
        _, flag = P.ccorr(*b, return_flag=True)
        rng = np.random.default_rng([SEED, 2, sum(map(ord, k))])
        nul_o[k] = P.swap_offset_null([b], [1.0], n_perm=n_perm, rng=rng)
        nul_r[k] = P.swap_null([b], [1.0], n_perm=n_perm, rng=rng)
        rows.append(dict(dataset=k, disease=TESTED[k], pairs=d.n, kind="counts" if d.counts else "array", clock_genes=len(d.wf.genes), gate_share=d.gate["share"],
                         gate_p=d.gate["p"], gated=bool(d.gate["p"] < 0.05 and d.n >= 20), coherence_flag=bool(flag), offset_null_sd_h=float(np.std(nul_o[k])),
                         offset_mde_h=MDE_K * float(np.std(nul_o[k])), r_null_sd=float(np.std(nul_r[k]))))
        say(k, "pairs", d.n, "gate p", round(d.gate["p"], 4))
    Pw = pd.DataFrame(rows)
    use = list(Pw[Pw.gated].dataset)
    w = np.array([D[k].n for k in use], float)
    pooled_null = P.swap_offset_null([B[k] for k in use], w, n_perm=n_perm, rng=np.random.default_rng([SEED, 3])) if use else np.array([np.nan])
    Pw["pooled_datasets"], Pw["pooled_offset_null_sd_h"], Pw["pooled_offset_mde_h"] = len(use), float(np.std(pooled_null)), MDE_K * float(np.std(pooled_null))
    Pw["H254.1_scored"] = bool(len(use) >= 3 and MDE_K * np.std(pooled_null) <= 1.0)
    Pw["own_test_scored"] = Pw.gated & (Pw.offset_mde_h <= 2.5)
    Pw.to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 250)
    print(Pw.round(4).to_string(), flush=True)
    if args.stage == "power":
        return

    comp = {kind: {k: [D[k].states(wf) for wf in D[k].pseudo(kind, n_sets, tag)] for k in D} for kind, tag in (("matched", 4), ("responsive", 5))}
    out, desc = [], []
    rng = np.random.default_rng([SEED, 6])

    def readout(k, d, b, nlo, nlr, cmp):
        c = P.ccorr(*b)
        off, r = float(np.angle(c) * HOURS), P.log_var_ratio(*b)
        n = len(b[0])
        bs = np.array([np.angle(P.ccorr(*[x[i] for x in b])) * HOURS for i in (rng.integers(0, n, n) for _ in range(n_boot))])
        bs = (bs - off + 12) % 24 - 12 + off
        row = dict(dataset=k, pairs=n, offset_h=off, lo90=float(np.nanquantile(bs, .05)), hi90=float(np.nanquantile(bs, .95)), coherence=float(abs(c)),
                   p_offset_two_sided=p_two(nlo, off), p_offset_delay=p_low(nlo, off), r=r, r_se=P.boot_se(*b, n_boot=n_boot, rng=rng), p_r_two_sided=p_two(nlr - np.nanmean(nlr), r - np.nanmean(nlr)))
        for kind, sets in cmp.items():
            po = np.array([np.angle(P.ccorr(*s)) * HOURS for s in sets])
            pr = np.array([P.log_var_ratio(*s) for s in sets])
            row.update({f"{kind}_abs_offset_q95": float(np.nanquantile(np.abs(po), .95)), f"{kind}_offset_q05": float(np.nanquantile(po, .05)),
                        f"{kind}_r_q05": float(np.nanquantile(pr, .05)), f"{kind}_r_q95": float(np.nanquantile(pr, .95)), f"{kind}_r_median": float(np.nanmedian(pr))})
        row["offset_clock_specific"] = bool(abs(off) > row["matched_abs_offset_q95"] and abs(off) > row["responsive_abs_offset_q95"])
        row["r_clock_specific"] = bool((r < row["matched_r_q05"] or r > row["matched_r_q95"]) and (r < row["responsive_r_q05"] or r > row["responsive_r_q95"]))
        b2 = d.states(drop=["PER1", "NR1D1"])
        row["offset_without_PER1_NR1D1_h"] = float(np.angle(P.ccorr(*b2)) * HOURS)
        return row

    for k, d in D.items():
        row = readout(k, d, B[k], nul_o[k], nul_r[k], {kind: comp[kind][k] for kind in comp})
        row.update(disease=TESTED[k], gated=bool(Pw.set_index("dataset").gated[k]), own_test_scored=bool(Pw.set_index("dataset").own_test_scored[k]))
        out.append(row)
    T = pd.DataFrame(out)
    m = T.own_test_scored
    T.loc[m, "p_offset_holm"] = holm(T.loc[m, "p_offset_two_sided"])
    T.loc[m, "p_r_holm"] = holm(T.loc[m, "p_r_two_sided"])
    pooled = {}
    if use:
        obs = P.pooled_offset([P.ccorr(*B[k]) for k in use], w)
        pooled = dict(datasets=";".join(use), pairs=int(w.sum()), pooled_offset_h=obs, null_sd_h=float(np.std(pooled_null)), p_delay=p_low(pooled_null, obs),
                      p_two_sided=p_two(pooled_null, obs), mde_h=MDE_K * float(np.std(pooled_null)), scored=bool(Pw["H254.1_scored"].iloc[0]))
        for kind in comp:
            po = np.array([P.pooled_offset([P.ccorr(*comp[kind][k][j]) for k in use], w) for j in range(n_sets)])
            pooled.update({f"{kind}_q05_h": float(np.nanquantile(po, .05)), f"{kind}_median_h": float(np.nanmedian(po)), f"{kind}_abs_q95_h": float(np.nanquantile(np.abs(po), .95))})
        pooled["clock_specific_delay"] = bool(obs < pooled["matched_q05_h"] and obs < pooled["responsive_q05_h"])
        t = T.set_index("dataset").loc[use]
        iv = 1 / t.r_se.to_numpy() ** 2
        pr = float(np.sum(iv * t.r) / iv.sum())
        nr = sum(wi * (nul_r[k] - np.nanmean(nul_r[k])) for wi, k in zip(iv / iv.sum(), use))
        pooled.update(pooled_r=pr, pooled_r_se=float(1 / np.sqrt(iv.sum())), p_r_two_sided=p_two(nr, pr - float(np.sum(iv / iv.sum() * [np.nanmean(nul_r[k]) for k in use]))))
        for kind in comp:
            v = np.array([np.nansum(iv / iv.sum() * [P.log_var_ratio(*comp[kind][k][j]) for k in use]) for j in range(n_sets)])
            pooled.update({f"{kind}_r_q05": float(np.quantile(v, .05)), f"{kind}_r_q95": float(np.quantile(v, .95))})
    for key, (name, ref_label) in DESC.items():
        d = Pair(key, ref_label, args.placebo)
        b = d.states()
        r2 = np.random.default_rng([SEED, 7, sum(map(ord, key))])
        row = readout(key, d, b, P.swap_offset_null([b], [1.0], n_perm=n_perm, rng=r2), P.swap_null([b], [1.0], n_perm=n_perm, rng=r2),
                      {kind: [d.states(wf) for wf in d.pseudo(kind, n_sets, tag)] for kind, tag in (("matched", 4), ("responsive", 5))})
        row.update(disease=name, gate_p=d.gate["p"])
        desc.append(row)
    T.to_csv(OUT / "score_datasets.csv", index=False)
    pd.DataFrame([pooled]).to_csv(OUT / "score_pooled.csv", index=False)
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    print(T.round(3).T.to_string(), pd.Series(pooled).to_string(), pd.DataFrame(desc).round(3).T.to_string(), sep="\n", flush=True)


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
