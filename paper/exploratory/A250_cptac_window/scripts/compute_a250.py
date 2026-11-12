"""A250: are a tumour's clock-controlled genes on the patient's body time while its drug-handling genes are not? (CPTAC-3)

Registered in paper/docs/registered_experiments.md (A250); confirmatory for the descriptive result of A249. Statistics and
gene sets are those of `compute_a249.py` (imported): S_T(set) = mean over genes of r(y of the tumour sample, z_N),
S_N(set) = mean of r(y of the normal sample, z_T), pooled over CPTAC types with patient-count weights; null: patients of
y permuted within type, one permutation for all genes, sets and sides. Delta(O, X) = S_T(O) - S_T(X).

Stages
  power  eligibility, null SD and MDE (2.49 x SD) of S_T, S_N and of each Delta, permuted patients only -> output/power.csv
  score  (after the power note is dated) everything registered -> output/score_{tests,contrasts,genes,pseudo,matched,descriptive}.csv
  `--placebo`: the patients of y are shuffled once within type before anything else (code test; output/placebo/).

    python compute_a250.py --stage power|score [--smoke] [--placebo]

Choices not fixed by the registration
- Blocks: the 7 types of the A186 pool; states from the full count tables with `clocksync.paired.cohort_states`, asserted
  equal to `a186_lib.states`. Duplicated gene symbols: the first column is kept.
- A gene's pooled r runs over the types where it is eligible; a set is scored with >= half of its genes (>= 100 for F)
  eligible in >= 4 types, and S averages the genes eligible in >= 1 type. Matched genes and pseudo-clock sets as in A249
  (pseudo sets drawn per type from the normals, set i pooled across types).
- Delta: one-sided p = share of permutations with Delta_null >= Delta; the bootstrap resamples patients within type with
  the same rows for both sets. Its MDE is 2.49 x the SD of Delta_null.
- H250.4: per-gene pooled r in the tumour sample against z_N, null from the same permutations; Holm over the three genes.
- H250.5: P249.0 and H249.1-2 with the thresholds of A249 (D = S_N - S_T of set C, 90 % bootstrap interval).
- Power rule: not scored if the MDE of S_T(O) > 0.042 or the MDE of Delta(O, C) > 0.045.
- Descriptive strata: G2M tertiles of A186 (`a186_lib.labels`), low and high, y and z re-centred within the stratum.
  No purity table of CPTAC is on disk: the purity strata are omitted.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXPL = HERE.parents[1]
sys.path.insert(0, str(EXPL / "A249_chrono_window" / "scripts"))
sys.path.insert(0, str(EXPL / "A186_cptac_c1" / "scripts"))
import compute_a249 as A  # noqa: E402
import a186_lib as L  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.paired import cohort_states, pseudo_waveforms  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
SEED = 250
POOL = ["LUAD", "CCRCC", "LSCC", "SKCM", "HNSCC", "PDAC", "UCEC"]
A.MIN_BLOCKS = 4
MDE_O, MDE_DELTA = 0.042, 0.045
GENES3 = ["NFIL3", "RORC", "ARNTL2"]
say = A.say


def cptac_block(tp, i, n_pseudo):
    cases, ct, cn = L.load_type(tp)
    keep = ~ct.columns.duplicated()
    ct, cn = ct.loc[:, keep], cn.loc[:, keep]
    ft, fn = ct.assign(depth=ct.sum(1)), cn.assign(depth=cn.sum(1))
    bt, _, bn, _ = cohort_states(ft, fn, True, A.WF)
    zT, zN = S.complex_state(bt), S.complex_state(bn)
    ref = L.states(tp, L.TYPES.index(tp), n_gate=20)
    assert np.allclose(zT - zT.mean(), ref["za"], atol=1e-6) and np.allclose(zN - zN.mean(), ref["zb"], atol=1e-6), tp
    ps = []
    if n_pseudo:
        lm = np.log(cn.div(cn.sum(1), axis=0).mean() + 1e-12)[(cn.mean() >= 1) | cn.columns.isin(list(A.WF.genes))]
        ps = [cohort_states(ft, fn, True, w) for w in pseudo_waveforms(lm, n_pseudo, np.random.default_rng([SEED, 1, i]), A.WF)]
    b = A.block(tp, A.lcpm(cn, fn.depth), A.lcpm(ct, ft.depth), cn, ct, zN, zT,
                [S.complex_state(p[2]) for p in ps], [S.complex_state(p[0]) for p in ps])
    say(tp, b["n"], "pairs", len(b["genes"]), "genes")
    return b


def subset(b, m):
    c = lambda v: v[m] - v[m].mean(0)  # noqa: E731
    return dict(b, n=int(m.sum()), stratum=np.zeros(int(m.sum()), int), yN=c(b["yN"]), yT=c(b["yT"]), zN=c(b["zN"]), zT=c(b["zT"]))


def boot_rows(blocks, rng):
    return [rng.choice(b["n"], b["n"]) for b in blocks]


def p_up(nul, obs):
    return float((1 + np.sum(nul >= obs)) / (1 + len(nul)))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o])
    out = np.empty(len(p))
    out[o] = np.minimum(adj, 1)
    return out


def run(args):
    sets = dict(A.SETS, F=A.fda_targets())
    excl = set(A.WF.genes) | {g for v in sets.values() for g in v}
    n_perm, n_pseudo, n_boot = (200, 5, 50) if args.smoke else (A.N_PERM, A.N_PSEUDO, A.N_BOOT)
    score = args.stage == "score"
    pool = POOL[:4] if args.smoke else POOL
    blocks = [cptac_block(tp, i, n_pseudo if score else 0) for i, tp in enumerate(pool)]
    if args.placebo:
        rng = np.random.default_rng([SEED, 8])
        for b in blocks:
            p = rng.permutation(b["n"])
            b["yN"], b["yT"] = b["yN"][p], b["yT"][p]
    perms = A.permutations(blocks, n_perm, [SEED, 3])
    n_tot = sum(b["n"] for b in blocks)
    nul = {(s, side): A.chunked_null(blocks, g, side, perms) for s, g in sets.items() for side in "NT"}
    say("nulls done")
    power = []
    for s, g in sets.items():
        e4, e1 = A.eligible(blocks, g)
        power.append(dict(item=f"set {s}", pairs=n_tot, types=len(blocks), genes=len(g), genes_in_4_types=e4, genes_scored=e1,
                          enters=bool(e4 >= (100 if s == "F" else len(g) / 2)),
                          **{f"{k}_S_{side}": f(nul[(s, side)]) for side in "TN"
                             for k, f in (("null_mean", np.mean), ("null_sd", np.std), ("mde", lambda v: A.MDE_K * np.std(v)))}))
    dnul = {x: nul[("O", "T")] - nul[(x, "T")] for x in "CFK"}
    for x in "CFK":
        power.append(dict(item=f"Delta(O,{x})", pairs=n_tot, types=len(blocks), null_mean_S_T=float(dnul[x].mean()),
                          null_sd_S_T=float(dnul[x].std()), mde_S_T=A.MDE_K * float(dnul[x].std())))
    P = pd.DataFrame(power)
    mde_o = float(P.loc[P.item == "set O", "mde_S_T"].iloc[0])
    mde_d = float(P.loc[P.item == "Delta(O,C)", "mde_S_T"].iloc[0])
    P["scored"] = bool(mde_o <= MDE_O and mde_d <= MDE_DELTA and P.enters.dropna().all())
    P.to_csv(OUT / "power.csv", index=False)
    print(P.round(4).to_string(), flush=True)
    if not score:
        return

    tests, prows, mrows, obs = [], [], [], {}
    for k, (s, g) in enumerate(sets.items()):
        ms = A.matched_sets(blocks, g, excl, np.random.default_rng([SEED, 4, k]))
        row = dict(set=s, pairs=n_tot, types=len(blocks), genes_scored=A.eligible(blocks, g)[1])
        for side in "TN":
            v, nl = float(A.set_stat(blocks, g, side)), nul[(s, side)]
            pse = A.set_stat(blocks, g, side, state_key="ZTp" if side == "N" else "ZNp")
            mat = np.array([A.set_stat(blocks, m, side) for m in ms])
            obs[(s, side)] = v
            row.update({f"S_{side}": v, f"null_mean_{side}": float(nl.mean()), f"null_sd_{side}": float(nl.std()), f"p_{side}": p_up(nl, v),
                        f"pseudo_median_{side}": float(np.median(pse)), f"pseudo_q95_{side}": float(np.quantile(pse, .95)),
                        f"p_clock_{side}": p_up(pse, v), f"matched_median_{side}": float(np.median(mat)),
                        f"matched_q95_{side}": float(np.quantile(mat, .95)), f"p_matched_{side}": p_up(mat, v),
                        f"mde_{side}": A.MDE_K * float(nl.std())})
            prows += [dict(set=s, side=side, i=i, value=float(x)) for i, x in enumerate(pse)]
            mrows += [dict(set=s, side=side, i=i, value=float(x)) for i, x in enumerate(mat)]
        tests.append(row)
        say("set", s, "scored")
    rng = np.random.default_rng([SEED, 7])
    boot = {k: [] for k in ("C", "F", "K", "D")}
    for i in range(n_boot):
        rows = boot_rows(blocks, rng)
        st = {s: float(A.set_stat(blocks, sets[s], "T", rows=rows)) for s in "OCFK"}
        for x in "CFK":
            boot[x].append(st["O"] - st[x])
        boot["D"].append(float(A.set_stat(blocks, sets["C"], "N", rows=rows)) - st["C"])
        if i % 500 == 0:
            say("bootstrap", i)
    con = []
    for x in "CFK":
        d, b = obs[("O", "T")] - obs[(x, "T")], np.array(boot[x])
        con.append(dict(contrast=f"Delta(O,{x})", value=d, null_mean=float(dnul[x].mean()), null_sd=float(dnul[x].std()),
                        p=p_up(dnul[x], d), lo90=float(np.quantile(b, .05)), hi90=float(np.quantile(b, .95)),
                        mde=A.MDE_K * float(dnul[x].std())))
    b = np.array(boot["D"])
    con.append(dict(contrast="D = S_N(C) - S_T(C)", value=obs[("C", "N")] - obs[("C", "T")], lo90=float(np.quantile(b, .05)),
                    hi90=float(np.quantile(b, .95))))
    C = pd.DataFrame(con)
    C.loc[C.contrast.isin(["Delta(O,F)", "Delta(O,K)"]), "p_holm"] = holm(C.loc[C.contrast.isin(["Delta(O,F)", "Delta(O,K)"]), "p"])

    grow, desc = [], []
    for s, g in sets.items():
        idx = A.index(blocks, g)
        nb = np.sum([ix >= 0 for ix in idx], axis=0)
        for side in "TN":
            zk = "zT" if side == "N" else "zN"
            per = [A.r_obs(b["y" + side], b[zk]) for b in blocks]
            cr = [A.unit(b["y" + side]).T @ np.conj(b[zk]) / np.sqrt(np.sum(abs(b[zk]) ** 2)) for b in blocks]
            pr = A.pooled(blocks, per, idx)
            pc = A.pooled(blocks, [np.real(c) for c in cr], idx) + 1j * A.pooled(blocks, [np.imag(c) for c in cr], idx)
            grow += [dict(set=s, side=side, gene=a, types=int(k), r=float(r),
                          implied_peak_h=float((-np.angle(c)) % (2 * np.pi) * S.HOURS) if k else np.nan)
                     for a, r, c, k in zip(g, pr, pc, nb)]
        for label, yk, zk in (("same_sample_normal", "yN", "zN"), ("same_sample_tumour", "yT", "zT")):
            desc.append(dict(item=label, set=s, value=float(np.nanmean(A.pooled(blocks, [A.r_obs(b[yk], b[zk]) for b in blocks], idx)))))
        for b in blocks:
            for side in "TN":
                desc.append(dict(item=f"S_{side}", set=s, type=b["name"], pairs=b["n"], value=float(A.set_stat([b], g, side))))
    G = pd.DataFrame(grow)
    g3 = []
    for a in GENES3:
        v = float(A.set_stat(blocks, [a], "T"))
        nl = A.chunked_null(blocks, [a], "T", perms)
        g3.append(dict(gene=a, r=v, null_mean=float(nl.mean()), null_sd=float(nl.std()), p=p_up(nl, v)))
    G3 = pd.DataFrame(g3)
    G3["p_holm"] = holm(G3.p)
    if not args.smoke:
        stad = cptac_block("STAD", 7, 0)
        for s in "OC":
            for side in "TN":
                desc.append(dict(item=f"S_{side}", set=s, type="STAD (not pooled)", pairs=stad["n"], value=float(A.set_stat([stad], sets[s], side))))
    for name, lab in (("G2M low", 0), ("G2M high", 2)):
        sub = [subset(b, L.labels(b["name"])[1] == lab) for b in blocks]
        for s in "OC":
            desc.append(dict(item=f"S_T, {name} tertile", set=s, pairs=sum(x["n"] for x in sub), value=float(A.set_stat(sub, sets[s], "T"))))
    T = pd.DataFrame(tests)
    T.to_csv(OUT / "score_tests.csv", index=False)
    C.to_csv(OUT / "score_contrasts.csv", index=False)
    G.to_csv(OUT / "score_genes.csv", index=False)
    G3.to_csv(OUT / "score_genes3.csv", index=False)
    pd.DataFrame(prows).to_csv(OUT / "score_pseudo.csv", index=False)
    pd.DataFrame(mrows).to_csv(OUT / "score_matched.csv", index=False)
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    print(T.round(4).T.to_string(), C.round(4).to_string(), G3.round(4).to_string(), sep="\n", flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--placebo", action="store_true")
    ap.add_argument("--smoke", action="store_true", help="4 types, 200 permutations, 5 pseudo sets, 50 bootstraps")
    args = ap.parse_args(argv)
    cli.setup(args)
    global OUT
    OUT = OUT / ("placebo" if args.placebo else "smoke" if args.smoke else "")
    OUT.mkdir(parents=True, exist_ok=True)
    run(args)


if __name__ == "__main__":
    main()
