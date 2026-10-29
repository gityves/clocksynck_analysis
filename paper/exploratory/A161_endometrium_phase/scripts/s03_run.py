"""A161: endometrial clock coherence, secretory vs proliferative (H161.1), pseudo-clock specificity, MDE and sensitivities (i)-(iv).
python -I s03_run.py [--nperm 10000] [--npseudo 200]
Coherence = clocksync.decompose.components (A87/A97), y = log2 expression z-scored per gene, centred within the sample group, no Poisson term.
T = coherence(secretory) - coherence(proliferative); H161.1: T < 0, labels permuted within endometriosis status, one-sided add-one p."""
import argparse, sys, re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
import a161_lib as L
import clocksync as cs
from clocksync._core.controls import CLOCK_LIKE
from clocksync.decompose import components
from clocksync.genes import canonical
from clocksync.waveform import Waveform

ap = argparse.ArgumentParser(); ap.add_argument("--nperm", type=int, default=10000); ap.add_argument("--npseudo", type=int, default=200)
args = ap.parse_args()
OUT, DER = L.OUT, L.EN / "derived"
WF = cs.gtex(); GENES = list(WF.genes)
rng = np.random.default_rng(161)
EPI, STR = ["EPCAM", "KRT8", "KRT18", "CDH1"], ["VIM", "PDGFRB", "COL1A1", "DCN"]
WARN = []


def wf_for(genes) -> Waveform:
    return Waveform(tuple(genes), WF.amp, WF.psi, None, None, "set")


def coh(Y, idx, wf) -> float:
    y = Y[idx] - Y[idx].mean(axis=0)
    return components(y, np.zeros(Y.shape[1]), 1, wf)["coherence"]


def Tstat(Y, wf, lab, idx=None, a=1, b=0):
    """coherence(label a) - coherence(label b) over the samples `idx`."""
    idx = np.arange(len(lab)) if idx is None else idx
    return coh(Y, idx[lab[idx] == a], wf) - coh(Y, idx[lab[idx] == b], wf)


def permuted(lab, strata, rng):
    out = lab.copy()
    for s in np.unique(strata):
        i = np.flatnonzero(strata == s); out[i] = lab[rng.permutation(i)]
    return out


def test(name, Yc, pseudo, lab, strata, nperm, a=1, b=0, mask=None):
    """lab: int labels over all rows of Yc; mask: rows analysed (default all with lab in {a,b}). Returns dict."""
    mask = np.isin(lab, [a, b]) if mask is None else mask
    rows = np.flatnonzero(mask); l = lab[rows]; st = np.asarray(strata)[rows]
    Y = Yc[rows]
    obs = Tstat(Y, WF, l, None, a, b)
    null = np.array([Tstat(Y, WF, permuted(l, st, rng), None, a, b) for _ in range(nperm)])
    r = dict(test=name, n_a=int((l == a).sum()), n_b=int((l == b).sum()), T=obs, coh_a=coh(Y, np.flatnonzero(l == a), WF),
             coh_b=coh(Y, np.flatnonzero(l == b), WF), p_less=cs.stats.perm_p(null, obs, "less", tie=0.0),
             null_mean=null.mean(), null_sd=null.std(ddof=1), MDE=2.49 * null.std(ddof=1), n_perm=nperm)
    if pseudo is not None:
        tp = np.array([Tstat(Yp[rows], w, l, None, a, b) for Yp, w in pseudo])
        r.update(pseudo_mean=tp.mean(), pseudo_q05=np.quantile(tp, .05), pseudo_median=np.median(tp),
                 pseudo_percentile=100 * np.mean(tp <= obs), clock_below_q05=bool(obs < np.quantile(tp, .05)), n_pseudo=len(tp))
        pd.DataFrame({"T_pseudo": tp}).to_csv(OUT / f"pseudoT_{re.sub('[^A-Za-z0-9]+', '_', name)}.csv", index=False)
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    return r


def zs(G: pd.DataFrame) -> pd.DataFrame:  # genes x samples -> z per gene over samples
    return G.sub(G.mean(axis=1), axis=0).div(G.std(axis=1, ddof=1), axis=0)


def residualise(Yfull, score_mat):
    X = np.column_stack([np.ones(len(score_mat)), score_mat])
    beta = np.linalg.lstsq(X, Yfull, rcond=None)[0]
    return Yfull - X @ beta


def draw_pseudo(Z: pd.DataFrame, npseudo, rng):
    """A127 variance-matched draw: pool = genes not clock-like (CLOCK_LIKE), not clock genes, not composition markers; for each clock
    gene the 20 nearest pool genes in standardised (mean, sd) of the raw log2 values (mean, sd given in Z.attrs); one drawn per gene."""
    G = Z.attrs["raw"]
    pool = [g for g in G.index if canonical(g) not in CLOCK_LIKE and g not in GENES and g not in EPI + STR]
    pm, ps = G.loc[pool].mean(axis=1), G.loc[pool].std(axis=1)
    z = lambda a, ref: (a - ref.mean()) / ref.std()
    P = np.column_stack([z(pm.values, pm), z(ps.values, ps)])
    tgt = np.column_stack([z(G.loc[GENES].mean(axis=1).values, pm), z(G.loc[GENES].std(axis=1).values, ps)])
    sets = []
    for _ in range(npseudo):
        used, s = set(), []
        for t in rng.permutation(len(GENES)):
            near = [j for j in np.argsort(((P - tgt[t]) ** 2).sum(1))[:40] if j not in used][:20]
            j = int(rng.choice(near)); used.add(j); s.append((t, pool[j]))
        sets.append([g for _, g in sorted(s)])
    mq = pd.DataFrame({"clock_mean": G.loc[GENES].mean(axis=1).values, "clock_sd": G.loc[GENES].std(axis=1).values,
                       "pseudo_mean_avg": [G.loc[[s[i] for s in sets]].mean(axis=1).mean() for i in range(14)],
                       "pseudo_sd_avg": [G.loc[[s[i] for s in sets]].std(axis=1).mean() for i in range(14)]}, index=GENES)
    return sets, mq


# =============================================== primary series GSE51981 ====================================================
meta = pd.read_csv(OUT / "gse51981_samples.csv").set_index("gsm")
G = pd.read_parquet(DER / "GSE51981_gene_log2.parquet")[meta.index]
assert all(g in G.index for g in GENES), "clock genes missing"
Z = zs(G); Z.attrs["raw"] = G
sets, mq = draw_pseudo(Z, args.npseudo, rng)
pd.DataFrame(sets).to_csv(OUT / "pseudo_sets_GSE51981.csv", index=False); mq.to_csv(OUT / "pseudo_match_quality.csv"); print(mq.round(2))
Yc = Z.loc[GENES].T.to_numpy()
pseudo = [(Z.loc[s].T.to_numpy(), wf_for(s)) for s in sets]
lab = meta.phase.map({"P": 0, "S": 1}).fillna(-1).astype(int).to_numpy()
endo = meta.endo.to_numpy()
res = []

# H161.1
res.append(test("H161.1 primary (S vs P, strata endo status)", Yc, pseudo, lab, endo, args.nperm))

# (i) composition
epi = Z.loc[[g for g in EPI if g in Z.index]].mean(axis=0).to_numpy(); stro = Z.loc[[g for g in STR if g in Z.index]].mean(axis=0).to_numpy()
print("composition genes found", [g for g in EPI + STR if g in Z.index])
S2 = np.column_stack([epi, stro])
Yc_r = residualise(Yc, S2)
pseudo_r = [(residualise(Yp, S2), w) for Yp, w in pseudo]
res.append(test("(i) composition-residualised", Yc_r, pseudo_r, lab, endo, args.nperm))
pd.DataFrame({"gsm": meta.index, "phase": meta.phase, "epi_score": epi, "stroma_score": stro}).to_csv(OUT / "composition_scores.csv", index=False)
cor = pd.DataFrame({"phase_S": (lab == 1).astype(float), "epi": epi, "stroma": stro})[lab >= 0].corr().round(3); cor.to_csv(OUT / "composition_vs_phase_corr.csv")

# (ii) controls without endometriosis or pelvic pathology
ctrl = (meta.source_name_ch1 == "Non-Endometriosis_No Uterine Pelvic Pathology").to_numpy()
res.append(test("(ii) controls without endometriosis/pathology", Yc, pseudo, lab, np.zeros(len(lab), int), args.nperm, mask=ctrl & (lab >= 0)))

# (iii) early / mid secretory vs proliferative
for sub in ("early_secretory", "mid_secretory"):
    l3 = np.where(meta.phase_label == "proliferative", 0, np.where(meta.phase_label == sub, 1, -1))
    res.append(test(f"(iii) {sub} vs proliferative (strata endo)", Yc, pseudo, l3, endo, args.nperm))

# (iv) replication: GSE6364 (GPL570, z within series), strata endometriosis / normal
m6 = pd.read_csv(OUT / "gse6364_samples.csv")
h6 = L.header(L.EN / "GSE6364_series_matrix.txt.gz")
h6["phase"] = h6["char"].str.extract(r"^(Proliferative|Early Secretory|Mid Secretory)")[0].map({"Proliferative": 0, "Early Secretory": 1, "Mid Secretory": 1})
h6["endo"] = h6["char"].str.contains("Endometriosis").astype(int)
G6 = pd.read_parquet(DER / "GSE6364_gene_log2.parquet")[h6.geo_accession]; Z6 = zs(G6)
print("GSE6364 counts", h6.phase.value_counts().to_dict())
res.append(test("(iv) GSE6364 (S vs P, strata endo status)", Z6.loc[GENES].T.to_numpy(), None, h6.phase.fillna(-1).astype(int).to_numpy(), h6.endo.to_numpy(), args.nperm))

# (iv) GSE141549: endometrium samples, cycle phase proliferative / secretory (medication, menstruation, unknown dropped), platform stratum
ys, labs, plats, rows = [], [], [], []
for tag, pf in (("HT12", "GPL10558"), ("WG6", "GPL13376")):
    h = L.header(L.EN / f"GSE141549-{pf}_series_matrix.txt.gz")
    h = h[h.tissue == "Endometrium"].copy()
    Gp = pd.read_parquet(DER / f"GSE141549_{tag}_gene_log2.parquet")
    h["col"] = h.title.str.strip()
    miss = [t for t in h.col if t not in Gp.columns]
    if miss: WARN.append(f"{tag}: {len(miss)} endometrium titles missing in xlsx: {miss[:5]}")
    h = h[h.col.isin(Gp.columns)]
    Gp = Gp[h.col]
    if not all(g in Gp.index for g in GENES):
        WARN.append(f"{tag}: clock genes missing {[g for g in GENES if g not in Gp.index]}"); continue
    Zp = zs(Gp)
    lp = h["cycle phase"].map({"proliferative": 0, "secretory": 1}).fillna(-1).astype(int).to_numpy()
    ys.append(Zp.loc[GENES].T.to_numpy()); labs.append(lp); plats += [pf] * len(h)
    rows += [dict(gsm=a, platform=pf, phase=b, replicate="Replicate" in t, no_stage=pd.isna(st))
             for a, b, t, st in zip(h.geo_accession, h["cycle phase"], h.title, h["disease stage"])]
pd.DataFrame(rows).to_csv(OUT / "gse141549_endometrium_samples.csv", index=False)
Y14 = np.vstack(ys); l14 = np.concatenate(labs)
print("GSE141549 phase counts", pd.Series(l14).value_counts().to_dict(), pd.crosstab(np.array(plats), l14).to_dict())
meta14 = pd.DataFrame(rows)
rep, nost = meta14.replicate.to_numpy(), meta14.no_stage.to_numpy()
# A: technical replicates dropped (same donors measured twice), all others kept
res.append(test("(iv) GSE141549 replicates dropped (S vs P, strata platform)", Y14, None, l14, np.array(plats), args.nperm, mask=np.isin(l14, [0, 1]) & ~rep))
# B: counts as registered (24 vs 43: replicates kept, the one secretory sample without disease stage dropped)
res.append(test("(iv) GSE141549 registered counts (S vs P, strata platform)", Y14, None, l14, np.array(plats), args.nperm, mask=np.isin(l14, [0, 1]) & ~nost))

out = pd.DataFrame(res); out.to_csv(OUT / "a161_results.csv", index=False)
(OUT / "a161_warnings.txt").write_text("\n".join(WARN))
print(out.T.to_string()); print(WARN)
