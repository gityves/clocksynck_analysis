"""A197 scoring (registered_experiments.md, A197, commit 6337b1b): does tumour clock coherence fall with ES-like dedifferentiation
beyond proliferation and subtype?
Within strata (breast: PAM50 x proliferation tertile; LUAD: proliferation x purity tertile; >= 30 tumours), T = mean over strata of
coh(high ES-like tertile) - coh(low); ES-like labels permuted within stratum. H197.1 SCAN-B, H197.2 METABRIC (Holm over two);
H197.3-4 GSE72094, GSE68465 (secondary, Holm over two). Coherence, tertiles, strata helpers and pseudo sets from A127 / A166.

Order: power (--power-only stops) -> tests, specificity, descriptives.
Run: python -I score_a197.py [--nperm 2000] [--npseudo 200] [--power-only]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A166_metabric" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import a166_lib as L  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

A = L.A  # a127_lib
OUT = HERE.parents[1] / "output"
GMT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/reference/c2.cgp.v2024.1.Hs.symbols.gmt")
K, LIMIT, SEED, MIN_STRATUM = 2.49, 0.15, 20261104, 30
TESTS = {"H197.1": "SCAN-B", "H197.2": "METABRIC", "H197.3": "GSE72094", "H197.4": "GSE68465"}


def say(*a):
    print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


def gmt(name):
    for line in GMT.read_text().splitlines():
        f = line.split("\t")
        if f[0] == name:
            return f[2:]
    raise KeyError(name)


EXCLUDE = set(L.sh.hallmark_g2m()) | set(gmt("BENPORATH_PROLIFERATION")) | set(gmt("BENPORATH_CYCLING_GENES")) | set(L.sh.PURITY) | set(L.GENES)
ES = [g for g in gmt("BENPORATH_ES_1") if g not in EXCLUDE and canonical(g) not in CLOCK_LIKE]
PRC2 = [g for g in gmt("BENPORATH_PRC2_TARGETS") if g not in EXCLUDE and canonical(g) not in CLOCK_LIKE]


def cohort(name):
    """(meta with pos, prolif, purity, pam50 (breast), es, prc2; expression samples x genes)."""
    if name == "SCAN-B":
        m, X = A.load()
    elif name == "METABRIC":
        m, X = L.load("metabric")
    else:
        d = EXP / "A171_luad_lihc" / "output"
        m, X = pd.read_parquet(d / f"{name.lower()}_meta.parquet"), pd.read_parquet(d / f"{name.lower()}_expr.parquet")
    m = m.copy()
    m["pos"] = np.arange(len(m))
    for k, gs in (("es", ES), ("prc2", PRC2)):
        g = [x for x in gs if x in X.columns]
        m[k] = L.sh._zmean(X, g).to_numpy()
        m[f"n_{k}_genes"] = len(g)
    return m, X


def strata(m, name, col="es", by_prolif=True):
    """(index arrays, `col` tertile labels within stratum)."""
    m = m.copy()
    if name in ("SCAN-B", "METABRIC"):
        m = m[m.pam50.notna()]
        m["ft"] = m.groupby("pam50").prolif.transform(L.tertile) if by_prolif else 0
        keys = ["pam50", "ft"]
    else:
        m["ft"] = L.tertile(m.prolif.to_numpy()) if by_prolif else 0
        m["pt"] = L.tertile(m.purity.to_numpy())
        keys = ["ft", "pt"]
    S, lab = [], []
    for _, g in m.groupby(keys, observed=True):
        if len(g) >= MIN_STRATUM:
            S.append(g.pos.to_numpy()); lab.append(L.tertile(g[col].to_numpy()))
    return L.Design("A197", "prolif", S, lab=lab)  # 'prolif' kind = high minus low tertile of the given labels


def Yof(X, genes):
    return L.LN2 * X[list(genes)].to_numpy(float)


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); n = len(p)
    out = np.empty(n); out[o] = np.maximum.accumulate((n - np.arange(n)) * p[o]).clip(max=1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
    ap.add_argument("--power-only", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    data, designs, power = {}, {}, []
    for j, (h, name) in enumerate(TESTS.items()):
        m, X = cohort(name)
        genes = [g for g in L.GENES if g in X.columns]
        Y, wf = Yof(X, genes), (A.WF if len(genes) == len(L.GENES) else A.WF.subset(genes))
        d = strata(m, name)
        null = d.null(Y, wf, a.nperm, np.random.default_rng([SEED, j]))
        data[h], designs[h] = (m, X, Y, wf, genes), d
        power.append(dict(test=h, cohort=name, tumours=d.n_tumours, strata=len(d.G), n_es_genes=int(m.n_es_genes.iloc[0]),
                          clock_genes=len(genes), null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                          mde=K * float(null.std(ddof=1)), limit=LIMIT, estimable=K * float(null.std(ddof=1)) <= LIMIT))
        data[h] += (null,)
        say(power[-1])
    pd.DataFrame(power).to_csv(OUT / "power.csv", index=False)
    if a.power_only:
        return

    rows, desc = [], []
    for j, (h, name) in enumerate(TESTS.items()):
        m, X, Y, wf, genes, null = data[h]
        d = designs[h]
        T = d.T(Y, wf)[0]
        p = float((1 + np.sum(null <= T)) / (1 + len(null)))
        pool = [g for g in L.pool_genes(X) if g not in set(ES) and g in X.columns]
        sets = L.pseudo_sets(X[pool + genes], pool, a.npseudo, np.random.default_rng([SEED, 10 + j])) if len(genes) == 14 else []
        pT = np.array([d.T(Yof(X, s), A.wf_for(s))[0] for s in sets]) if sets else np.array([np.nan])
        pw = next(r for r in power if r["test"] == h)
        rows.append(dict(test=h, cohort=name, T=T, p_one_sided=p, null_mean=pw["null_mean"], mde=pw["mde"], estimable=pw["estimable"],
                         pseudo_q05=float(np.nanquantile(pT, .05)), pseudo_median=float(np.nanmedian(pT)),
                         pseudo_percentile=float(np.nanmean(pT <= T) * 100), clock_specific=bool(T < np.nanquantile(pT, .05)),
                         tumours=d.n_tumours, strata=len(d.G)))
        # descriptives
        r = dict(cohort=name, rho_es_prolif=spearmanr(m.es, m.prolif)[0], rho_es_prc2=spearmanr(m.es, m.prc2)[0],
                 T_without_prolif_strata=strata(m, name, by_prolif=False).T(Y, wf)[0],
                 T_prc2_high_minus_low=strata(m, name, "prc2").T(Y, wf)[0])
        if "er" in m:
            r["rho_es_er"] = spearmanr(m.es, m.er, nan_policy="omit")[0]
        loo = [d.T(Yof(X, [g for g in genes if g != x]), A.WF.subset([g for g in genes if g != x]))[0] for x in genes]
        r.update(T_leave_one_gene_out_min=float(np.min(loo)), T_leave_one_gene_out_max=float(np.max(loo)))
        desc.append(r)
    t = pd.DataFrame(rows)
    t["p_holm"] = np.nan
    for fam in (["H197.1", "H197.2"], ["H197.3", "H197.4"]):
        k = t.test.isin(fam) & t.estimable
        if k.any():
            t.loc[k, "p_holm"] = holm(t.loc[k, "p_one_sided"])
    t["verdict"] = np.where(~t.estimable, "not estimable", np.where(t.p_holm < .05, "pass", "fail"))
    t.to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(desc).to_csv(OUT / "descriptive.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.round(4).to_string(index=False)); print(pd.DataFrame(desc).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
