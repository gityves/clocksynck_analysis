"""A268 (registered in `paper/docs/registered_experiments.md`): does a glucocorticoid or temperature reader keep a fixed phase
relation to the host clock state across tissues and people, clock-specifically, and is it confirmed in living people?

Readers (scout T4, relaxed rule chosen by the user): glucocorticoid FKBP5 + TSC22D3; cold CIRBP; heat HSP90AA1. Reader value
per sample = mean over its genes of the within-tissue z of log expression.
Host (same construction in both cohorts): per sample, the clock-plane state of the GTEx waveform (unweighted projector on
gene-centred log expression, `clocksync.paired.projector(counts=False)`), centred within tissue; the host of a (person,
tissue) is the mean state of the person's OTHER tissues (leave-one-tissue-out). Variant: waveform without PER1.
Statistic per reader: per tissue k the complex loading beta_k = sum y conj(z) / sum |z|^2 of the centred reader on the
centred host; consistency K = |sum_k beta_k| / sum_k |beta_k|; p against 10,000 permutations of the reader across persons
within tissue. Companions: 200 readers of expression-matched genes (same gene count, `pseudo_clock_genes` bins on the GTEx
mean log CPM, CLOCK_LIKE excluded) and 200 of responsive genes (`RESPONSIVE_GENES` minus the reader's genes, same count).
- GTEx: fast deaths, the 15 host tissues of the cue line (`gtex_cues/_shared.TISSUES`); body tissues = the 12 outside brain
  and pituitary carry K; brain (cerebellum, cortex) descriptive: projection of each brain beta on the unit vector of
  sum_body beta, relative to the median body projection, 200 donor bootstraps.
- MuTHER (living confirmation; fat and skin, the A255 matrices): beta_fat on the skin host and beta_skin on the fat host;
  pooled beta = their sum; p against 10,000 permutations of persons within tissue; argument difference to the GTEx body
  argument arg(sum_body beta).

    python score_a268.py [--jobs 6]
Outputs (`../output/`): tests.csv (one row per reader and host variant), companions.csv, brain.csv, muther.csv.
Controlled access: GTEx time of death is not used; outputs are aggregates.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_3/scripts/gtex_cues"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE, RESPONSIVE_GENES, pseudo_clock_genes  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, N_PERM, N_SETS, N_BOOT = 268, 10_000, 200, 200
READERS = {"glucocorticoid": ("FKBP5", "TSC22D3"), "cold": ("CIRBP",), "heat": ("HSP90AA1",)}
WF = cs.gtex()
WF_NOPER1 = WF.drop(["PER1"])
BRAIN = ("Brain - Cerebellum", "Brain - Cortex")
BODY = tuple(t for t in S.TISSUES if t not in BRAIN and t != "Pituitary")
PREP = PAPER / "exploratory/A254_paired_lesions/output/prep"
G = {}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def states(logx: pd.DataFrame, tissue, wf) -> np.ndarray:
    """Centred clock-plane state (complex) of the rows of one tissue (log expression, persons x genes)."""
    g = list(wf.genes)
    y = logx[g] - logx[g].mean()
    L = P.projector(y, wf, counts=False)
    b = y.to_numpy(float) @ L.T
    z = b[:, 0] + 1j * b[:, 1]
    return z - z.mean()


def hosts(tables: dict, wf) -> dict:
    """tissue -> Series (person -> leave-one-tissue-out mean state of the person's other tissues)."""
    st = {k: pd.Series(states(v, k, wf), index=v.index) for k, v in tables.items()}
    allz = pd.DataFrame(st)
    out = {}
    for k in tables:
        other = allz.drop(columns=k)
        out[k] = other.mean(axis=1, skipna=True)[other.notna().sum(1) >= 1]
    return out


def reader_values(tables: dict, genes) -> dict:
    out = {}
    for k, v in tables.items():
        x = v[list(genes)]
        out[k] = ((x - x.mean()) / x.std()).mean(axis=1)
    return out


def loadings(rd: dict, hs: dict, tissues, perm_rng=None) -> dict:
    b = {}
    for k in tissues:
        y, z = rd[k], hs[k]
        idx = y.index.intersection(z.index)
        yy, zz = y.loc[idx].to_numpy(float), z.loc[idx].to_numpy()
        yy, zz = yy - yy.mean(), zz - zz.mean()
        if perm_rng is not None:
            yy = yy[perm_rng.permutation(len(yy))]
        b[k] = np.sum(yy * np.conj(zz)) / np.sum(np.abs(zz) ** 2)
    return b


def K(b: dict) -> float:
    v = np.array(list(b.values()))
    return float(np.abs(v.sum()) / np.abs(v).sum())


def k_test(rd, hs, tissues, n_perm, rng):
    obs = K(loadings(rd, hs, tissues))
    null = np.array([K(loadings(rd, hs, tissues, rng)) for _ in range(n_perm)])
    return obs, float((1 + np.sum(null >= obs)) / (1 + n_perm))


def companion_K(args):
    genes, variant = args
    rd = reader_values(G["gtex"], genes)
    return K(loadings(rd, G["hosts"][variant], BODY))


def gtex_tables():
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(a.var_names)
    a.file.close()
    pb, _ = gtex.load(tissues=list(S.TISSUES), hardy=(1, 2), genes=genes, depth="library_size", require_tod=False)
    tables, lm = {}, None
    X = pb[genes].to_numpy(np.float32)
    lg = np.log(X / pb.depth.to_numpy(np.float32)[:, None] * 1e6 + 0.5)
    lm = pd.Series(np.log((X / pb.depth.to_numpy(np.float32)[:, None]).mean(0) + 1e-12), index=genes)[X.mean(0) >= 1]
    for k in S.TISSUES:
        m = (pb.compartment == k).to_numpy()
        tables[k] = pd.DataFrame(lg[m], index=pb.donor[m].to_numpy(), columns=genes)
    return tables, lm


def muther_power(n_draw=200, n_perm=1000, r2s=(0.0, 0.014, 0.03)):
    """Pass rate of the MuTHER confirmation test with an inserted reader-host relation (reader permuted within tissue,
    plus a relation of random fixed phase scaled to R^2 of the reader); per reader."""
    mu = (PREP / "muther_fat_expr.parquet", PREP / "muther_skin_expr.parquet")
    M = {k: pd.read_parquet(f).astype(float) * np.log(2.0) for k, f in zip(("fat", "skin"), mu)}
    M = {k: v.loc[:, ~v.columns.duplicated()] for k, v in M.items()}
    mh = hosts(M, WF.subset([g for g in WF.genes if g in M["fat"] and g in M["skin"]]))
    rows = []
    for r, genes in READERS.items():
        base = reader_values(M, genes)
        for r2 in r2s:
            hits = 0
            for i in range(n_draw):
                rng = np.random.default_rng([SEED, 5, i, int(r2 * 1000)])
                ph = rng.uniform(0, 2 * np.pi)
                rd = {}
                for k in ("fat", "skin"):
                    y, z = base[k], mh[k]
                    idx = y.index.intersection(z.index)
                    v = rng.permutation(y.loc[idx].to_numpy(float))
                    if r2 > 0:
                        sg = np.real(np.exp(-1j * ph) * (z.loc[idx].to_numpy() - z.loc[idx].mean()))
                        v = v + sg / sg.std() * v.std() * np.sqrt(r2 / (1 - r2))
                    rd[k] = pd.Series(v, index=idx)
                b = loadings(rd, mh, ("fat", "skin"))
                obs = abs(b["fat"] + b["skin"])
                null = [abs(sum(loadings(rd, mh, ("fat", "skin"), rng).values())) for _ in range(n_perm)]
                hits += (1 + np.sum(np.array(null) >= obs)) / (1 + n_perm) < 0.05
            rows.append(dict(reader=r, r2=r2, draws=n_draw, pass_rate=hits / n_draw))
            log(rows[-1])
    pd.DataFrame(rows).to_csv(OUT / "power_muther.csv", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--power", action="store_true", help="MuTHER power only (before scoring)")
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.power:
        return muther_power()
    rng = np.random.default_rng([SEED, 0])
    G["gtex"], lm = gtex_tables()
    log("GTEx tables", {k: len(v) for k, v in G["gtex"].items()})
    G["hosts"] = {"all": hosts(G["gtex"], WF), "noPER1": hosts(G["gtex"], WF_NOPER1)}
    mu = (PREP / "muther_fat_expr.parquet", PREP / "muther_skin_expr.parquet")
    M = {k: pd.read_parquet(f).astype(float) * np.log(2.0) for k, f in zip(("fat", "skin"), mu)}
    M = {k: v.loc[:, ~v.columns.duplicated()] for k, v in M.items()}
    MH = {"all": hosts(M, WF.subset([g for g in WF.genes if g in M["fat"] and g in M["skin"]])),
          "noPER1": hosts(M, WF_NOPER1.subset([g for g in WF_NOPER1.genes if g in M["fat"] and g in M["skin"]]))}
    rows, comp, brain, mrows = [], [], [], []
    for r, genes in READERS.items():
        rd = reader_values(G["gtex"], genes)
        for variant in ("all", "noPER1"):
            hs = G["hosts"][variant]
            k_obs, p = k_test(rd, hs, BODY, N_PERM, rng)
            b = loadings(rd, hs, BODY)
            body_arg = float(np.angle(sum(b.values())) * 12 / np.pi)
            # companions (gene sets of the reader's size)
            sets = {"matched": pseudo_clock_genes(lm, N_SETS, np.random.default_rng([SEED, 1, len(rows)]), genes=list(genes)),
                    "responsive": []}
            pool = [g for g in RESPONSIVE_GENES if g not in genes and g in lm.index]
            rr = np.random.default_rng([SEED, 2, len(rows)])
            sets["responsive"] = [list(rr.choice(pool, len(genes), replace=False)) for _ in range(N_SETS)]
            q = {}
            for kind, ss in sets.items():
                with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pl:
                    ks = np.array(pl.map(companion_K, [(s, variant) for s in ss]))
                comp += [dict(reader=r, host=variant, kind=kind, set=i, K=x) for i, x in enumerate(ks)]
                q[kind] = (float(np.quantile(ks, 0.95)), float(np.median(ks)), float(100 * np.mean(ks < k_obs)))
            # brain (descriptive)
            u = sum(b.values()) / abs(sum(b.values()))
            proj = {k: float(np.real(v * np.conj(u))) for k, v in loadings(rd, hs, S.TISSUES).items()}
            body_med = float(np.median([proj[k] for k in BODY]))
            bs = []
            donors = sorted(set().union(*[set(rd[k].index) for k in S.TISSUES]))
            for _ in range(N_BOOT):
                pick = rng.choice(donors, len(donors))
                hsb_aligned = {}
                for k in S.TISSUES:
                    y, z = rd[k].reindex(pick), hs[k].reindex(pick)
                    ok = (y.notna() & z.notna()).to_numpy()
                    hsb_aligned[k] = (y.to_numpy()[ok], z.to_numpy()[ok])
                bb = {k: np.sum((y - y.mean()) * np.conj(z - z.mean())) / np.sum(np.abs(z - z.mean()) ** 2)
                      for k, (y, z) in hsb_aligned.items()}
                ub = sum(bb[k] for k in BODY)
                ub = ub / abs(ub)
                pr = {k: float(np.real(v * np.conj(ub))) for k, v in bb.items()}
                med = np.median([pr[k] for k in BODY])
                bs.append([pr[k] / med for k in BRAIN])
            bs = np.array(bs)
            for i, k in enumerate(BRAIN):
                brain.append(dict(reader=r, host=variant, tissue=k, ratio_to_body=proj[k] / body_med,
                                  lo=float(np.quantile(bs[:, i], 0.025)), hi=float(np.quantile(bs[:, i], 0.975))))
            # MuTHER
            mh = MH[variant]
            mrd = reader_values(M, [g for g in genes if g in M["fat"]])
            mb = loadings(mrd, mh, ("fat", "skin"))
            pooled = mb["fat"] + mb["skin"]
            mnull = []
            for _ in range(N_PERM):
                x = loadings(mrd, mh, ("fat", "skin"), rng)
                mnull.append(abs(x["fat"] + x["skin"]))
            m_p = float((1 + np.sum(np.array(mnull) >= abs(pooled))) / (1 + N_PERM))
            arg_m = float(np.angle(pooled) * 12 / np.pi)
            diff = float((arg_m - body_arg + 12) % 24 - 12)
            mrows.append(dict(reader=r, host=variant, beta_fat_abs=abs(mb["fat"]), beta_skin_abs=abs(mb["skin"]),
                              pooled_abs=abs(pooled), pooled_arg_h=arg_m, p=m_p, gtex_body_arg_h=body_arg,
                              arg_diff_h=diff, confirmed=bool(m_p < 0.05 and abs(diff) <= 3)))
            row = dict(reader=r, genes="+".join(genes), host=variant, tissues=len(BODY), K=k_obs, p=p, body_arg_h=body_arg,
                       matched_q95=q["matched"][0], matched_median=q["matched"][1], matched_pct=q["matched"][2],
                       responsive_q95=q["responsive"][0], responsive_median=q["responsive"][1],
                       responsive_pct=q["responsive"][2])
            row["specific"] = bool(k_obs > q["matched"][0] and k_obs > q["responsive"][0])
            rows.append(row)
            log(r, variant, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})
    T = pd.DataFrame(rows)
    T["p_holm"] = np.nan
    for v in ("all", "noPER1"):
        m = T.host == v
        ps = T.loc[m, "p"].to_numpy()
        o = np.argsort(ps)
        adj = np.minimum(1, np.maximum.accumulate((len(ps) - np.arange(len(ps))) * ps[o]))
        h = np.empty_like(adj)
        h[o] = adj
        T.loc[m, "p_holm"] = h
    Mt = pd.DataFrame(mrows)
    T = T.merge(Mt[["reader", "host", "p", "arg_diff_h", "confirmed"]].rename(columns={"p": "muther_p"}), on=["reader", "host"])
    T["verdict"] = np.where((T.p_holm < 0.05) & T.specific & T.confirmed, "PASS (confirmed)",
                   np.where((T.p_holm < 0.05) & T.specific, "GTEx only (not confirmed)",
                   np.where(T.p_holm < 0.05, "not specific", "fail")))
    T.to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(comp).to_csv(OUT / "companions.csv", index=False)
    pd.DataFrame(brain).to_csv(OUT / "brain.csv", index=False)
    Mt.to_csv(OUT / "muther.csv", index=False)
    pd.set_option("display.width", 250)
    print(T.round(4).to_string(index=False), pd.DataFrame(brain).round(3).to_string(index=False),
          Mt.round(4).to_string(index=False), sep="\n", flush=True)


if __name__ == "__main__":
    main()
