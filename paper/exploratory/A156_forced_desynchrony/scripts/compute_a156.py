"""A156 (registered 2026-10-28): does the untimed core-cue offset detect forced desynchrony in human blood (GSE48113)?

Estimator (A124 / A115, arrays convention): log2 values, equal gene weights, one projector per gene set (core = GTEx 14 genes,
cue and locked = A81b sets) over the analysis arrays; per subject and day offset = arg of the complex correlation (`S.ccorr`,
states centred within subject x day) of the core state with the cue (locked) state. Delta_i = wrap(offset D4 - offset D1).
H156.1: |circular mean of Delta_i| (unweighted); H156.2: |mean d core-cue| - |mean d core-locked| (one-sided). Null: day labels
swapped within subject, which negates Delta_i exactly (per-subject-day offsets do not depend on labels): all 2^n sign
vectors if n <= 16, else 10,000 random ones (identity included in the count: p = (1 + #{null >= obs}) / (1 + N)).
Holm over the two. Seed 156.

    python compute_a156.py        (sccc python; reads output/samples.csv, probe_map.csv.gz and the series matrix)
Outputs: output/*.csv
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER.parent)); sys.path.insert(0, str(PAPER)); sys.path.insert(0, str(HERE))
sys.path.insert(0, str(PAPER / "exploratory/A115_A116_psoriasis/scripts"))
import _lib as L
import clocksync as cs
from clocksync import stats as S
from clocksync.config import limit_threads
from clocksync.genes import canonical
from _common import EXTRA

limit_threads(1)
OUT = HERE.parent / "output"
GEO = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/geo/GSE48113")
HOURS = S.HOURS
W = 2 * np.pi / 24
NB, NPERM, SEED = 2000, 10000, 156

smp = pd.read_csv(OUT / "samples.csv", dtype=str, keep_default_na=False)
smp["day"] = smp["day sample was taken in protocol"].astype(int)
smp["tp"] = smp["time point"].astype(int)
t = smp["time sample taken"].str.extract(r"^(\d+):(\d+)$").astype(float)
smp["clock_h"] = t[0] + t[1] / 60  # NaN for BB0264 point 7
smp["subject_day"] = smp.subject + "|D" + smp.day.astype(str)
pm = pd.read_csv(OUT / "probe_map.csv.gz", keep_default_na=False)
pm = pm[pm.control_type == "FALSE"]
expr = pd.read_csv(GEO / "GSE48113_series_matrix.txt.gz", sep="\t", comment="!", index_col=0, na_values=["", "null"])
expr.index = expr.index.astype(str)
expr = expr[smp.gsm.tolist()]

sets = {"core": cs.gtex()}
ref = pd.read_csv(PAPER / "results/figure_3/output/gtex_cues/a81b_sets.csv")
for k in ("cue", "locked"):
    x = ref[ref.set == k]
    sets[k] = cs.Waveform(tuple(x.gene), x.amp.to_numpy(), x.acro_h.to_numpy() * 2 * np.pi / 24, None, None, k)
ALL_GENES = sorted({g for w in sets.values() for g in w.genes})


def gene_table(alias: bool):
    """Array x gene table (mean of log2 over the gene's non-control probes present) and the mapping table."""
    canon = (lambda s: EXTRA.get(canonical(s), canonical(s))) if alias else canonical
    psym = pm.symbol.map(lambda s: canon(s) if s else "")
    cols, info = {}, []
    for g in ALL_GENES:
        pr = pm.probe[psym == canon(g)].astype(str)
        pr = pr[pr.isin(expr.index)]
        info.append(dict(gene=g, n_probes=len(pr)))
        if len(pr):
            cols[g] = expr.loc[pr].mean(axis=0, skipna=True).to_numpy()
    return pd.DataFrame(cols, index=smp.index), pd.DataFrame(info)


def impute(X, rows):
    """Gene missing in an array (all probes missing): mean of that gene within the subject-day (then within the day)."""
    X = X.loc[rows].copy(); n = int(X.isna().to_numpy().sum()); fb = 0
    for g in X.columns[X.isna().any()]:
        m = X[g].groupby(smp.loc[rows, "subject_day"]).transform("mean")
        X[g] = X[g].fillna(m)
        if X[g].isna().any():
            fb += int(X[g].isna().sum()); X[g] = X[g].fillna(X[g].groupby(smp.loc[rows, "day"]).transform("mean"))
    return X, n, fb


def circ_sd(x):
    return float(np.sqrt(-2 * np.log(abs(np.mean(np.exp(1j * np.asarray(x) / HOURS))))) * HOURS)


def cmean(x):
    return float(np.angle(np.mean(np.exp(1j * np.asarray(x) / HOURS))) * HOURS)


def cmean_abs(x):
    return abs(cmean(x))


def holm(p):
    o = np.argsort(p); out = np.empty(len(p)); run = 0
    for r, i in enumerate(o):
        run = max(run, (len(p) - r) * p[i]); out[i] = min(run, 1)
    return out


def analyse(label, subjects, alias, rng, full=False):
    X, mapinfo = gene_table(alias)
    rows = smp.index[smp.subject.isin(subjects)]
    Xi, n_imp, n_fb = impute(X, rows)
    s = smp.loc[rows]
    st, ng = {}, {}
    for k, w in sets.items():
        (b, _), ng[k] = L.states_for(Xi.assign(depth=1.0), w, counts=False)
        st[k] = b[:, 0] + 1j * b[:, 1]
    sd = s.subject_day.to_numpy()
    z = {}  # (subject, day) -> complex correlations
    for sub in subjects:
        for d in (1, 4):
            m = (s.subject == sub).to_numpy() & (s.day == d).to_numpy()
            if m.sum() < 3:
                raise ValueError(sub)
            b = {k: np.c_[st[k][m].real, st[k][m].imag] for k in st}
            z[sub, d] = {k: S.ccorr(b["core"], b[k]) for k in ("cue", "locked")}
    D = {k: np.array([S.wrap((np.angle(z[u, 4][k]) - np.angle(z[u, 1][k])) * HOURS) for u in subjects]) for k in ("cue", "locked")}
    n = len(subjects)
    obs1 = cmean_abs(D["cue"]); obs2 = obs1 - cmean_abs(D["locked"])
    if n <= 16:
        signs = np.array(list(itertools.product((1, -1), repeat=n))); exact = True
    else:
        signs = rng.choice([1, -1], size=(NPERM, n)); exact = False
    def mean_abs(M):  # circular-mean magnitude per sign vector (rows of M are sign vectors)
        return np.abs(np.angle((np.exp(1j * M / HOURS)).mean(1)) * HOURS)
    a1 = mean_abs(signs * D["cue"][None]); a2 = a1 - mean_abs(signs * D["locked"][None])
    if exact:
        p1, p2 = float(np.mean(a1 >= obs1 - 1e-12)), float(np.mean(a2 >= obs2 - 1e-12))
    else:
        p1 = (1 + float(np.sum(a1 >= obs1 - 1e-12))) / (1 + len(a1)); p2 = (1 + float(np.sum(a2 >= obs2 - 1e-12))) / (1 + len(a2))
    sgn = cmean(signs * 0 + 0) if False else None
    signed_null = np.angle(np.exp(1j * (signs * D["cue"][None]) / HOURS).mean(1)) * HOURS
    sd_null = circ_sd(signed_null)
    ph = holm([p1, p2])
    res = dict(analysis=label, n_subjects=n, n_arrays=len(rows), n_null=len(a1), exact=exact, genes_core=ng["core"],
               genes_cue=ng["cue"], genes_locked=ng["locked"], imputed_gene_arrays=n_imp, imputed_fallback=n_fb,
               d_core_cue_mean_h=cmean(D["cue"]), d_core_locked_mean_h=cmean(D["locked"]),
               H1_stat=obs1, H1_p=p1, H1_holm=ph[0], H2_stat=obs2, H2_p=p2, H2_holm=ph[1],
               null_circ_sd_h=sd_null, mde80_h=(1.96 + 0.8416) * sd_null, null_q95_abs_h=float(np.quantile(a1, 0.95)),
               p_min=1 / len(a1) if exact else 1 / (1 + len(a1)))
    return (res, dict(D=D, z=z, st=st, s=s, mapinfo=mapinfo, rows=rows, Xi=Xi, signed_null=signed_null)) if full else res


cnt = smp.groupby(["subject", "day"]).size().unstack(fill_value=0)
prim = sorted(cnt.index[(cnt[1] >= 5) & (cnt[4] >= 5)])
s77 = sorted(cnt.index[(cnt[1] == 7) & (cnt[4] == 7)])
s21 = sorted(cnt.index[(cnt[1] >= 1) & (cnt[4] >= 1)])
rng = np.random.default_rng(SEED)
res, aux = analyse("primary (>=5 per day, canonical+8 aliases)", prim, True, rng, full=True)
pd.DataFrame([res]).T.rename(columns={0: "value"}).to_csv(OUT / "primary_summary.csv")
tests = pd.DataFrame([
    dict(test="H156.1", statistic="|circ mean of d core-cue| (h)", observed=res["H1_stat"], p=res["H1_p"], p_holm=res["H1_holm"], alt="two-sided"),
    dict(test="H156.2", statistic="|d core-cue| - |d core-locked| (h)", observed=res["H2_stat"], p=res["H2_p"], p_holm=res["H2_holm"], alt="greater")])
tests["pass"] = tests.p_holm < 0.05
tests.to_csv(OUT / "tests.csv", index=False)
pd.DataFrame([dict(null_circ_sd_h=res["null_circ_sd_h"], mde80_h=res["mde80_h"], null_q95_abs_h=res["null_q95_abs_h"],
                   n_null=res["n_null"])]).to_csv(OUT / "mde.csv", index=False)

# sensitivities (descriptive companions; same code)
sens = [analyse("15 subjects 7+7 (exact), canonical+8 aliases", s77, True, rng),
        analyse("all 21 subjects, canonical+8 aliases", s21, True, rng),
        analyse("primary subjects, canonical only", prim, False, rng)]
pd.DataFrame(sens).to_csv(OUT / "sensitivities.csv", index=False)
aux["mapinfo"].to_csv(OUT / "gene_probe_counts_alias.csv", index=False)
_X, _ = gene_table(True)
_X = _X.loc[aux["rows"]]
_miss = _X.isna().groupby(smp.loc[aux["rows"], "subject_day"]).all()
pd.DataFrame(dict(gene=_X.columns, arrays_missing=_X.isna().sum().to_numpy(), subject_days_all_missing=_miss.sum().to_numpy())
             ).query("arrays_missing > 0").to_csv(OUT / "imputation_by_gene_primary.csv", index=False)
cov = []
for lab, al in (("canonical+8 aliases", True), ("canonical only", False)):
    _, mi = gene_table(al)
    for k, w in sets.items():
        cov.append(dict(mapping=lab, set=k, genes=len(w.genes), mapped=int(mi[mi.gene.isin(w.genes)].n_probes.gt(0).sum())))
pd.DataFrame(cov).to_csv(OUT / "gene_coverage.csv", index=False)

# ---- descriptives (primary set) ---------------------------------------------------------------------------------
D, z, st, s = aux["D"], aux["z"], aux["st"], aux["s"]
subs = prim
rb = np.random.default_rng(SEED + 1)
pd.DataFrame({"subject": subs, "d_core_cue_h": D["cue"], "d_core_locked_h": D["locked"]}).assign(
    subject=[f"S{i + 1:02d}" for i in range(len(subs))]).to_csv(OUT / "per_subject_delta_anonymised.csv", index=False)


def ci(fn, v):
    v = np.asarray(v); b = [fn(v[rb.integers(0, len(v), len(v))]) for _ in range(NB)]
    c = fn(v); w = S.wrap(np.asarray(b) - c)
    return c, float(c + np.quantile(w, 0.025)), float(c + np.quantile(w, 0.975))


# coherence per day
coh = []
for d in (1, 4):
    for k in ("cue", "locked"):
        zz = np.array([z[u, d][k] for u in subs])
        coh.append(dict(day=d, pair=f"core-{k}", mean_abs_z=float(np.abs(zz).mean()), abs_mean_z=float(abs(zz.mean())),
                        circ_mean_offset_h=cmean(np.angle(zz) * HOURS)))
pd.DataFrame(coh).to_csv(OUT / "coherence_per_day.csv", index=False)

# deltas with bootstrap
rows_d = []
for k in ("cue", "locked"):
    c, lo, hi = ci(cmean, D[k]); rows_d.append(dict(item=f"mean d core-{k} (h)", estimate=c, lo=lo, hi=hi, n=len(subs)))
pd.DataFrame(rows_d).to_csv(OUT / "descriptive_deltas.csv", index=False)

# timed cosinor: state ~ A exp(i (w t + theta)); theta (h) = arg sum (c_t - mean) exp(-i w t); needs >= 4 timed samples
th = {}
for sub in subs:
    for d in (1, 4):
        m = ((s.subject == sub) & (s.day == d) & s.clock_h.notna()).to_numpy()
        e = np.exp(-1j * W * s.clock_h.to_numpy()[m])
        rec = {"n": int(m.sum())}
        for k in ("core", "cue", "locked"):
            c = st[k][m]; v = ((c - c.mean()) * e).sum()
            rec[k] = float(np.angle(v) * HOURS); rec[k + "_amp"] = float(abs(v) / m.sum())
        th[sub, d] = rec
tim = pd.DataFrame([dict(subject=u, day=d, **r) for (u, d), r in th.items()])
tim["subject"] = tim.subject.map({u: f"S{i + 1:02d}" for i, u in enumerate(subs)})
tim.to_csv(OUT / "cosinor_per_subject_day_anonymised.csv", index=False)
dphi = {k: np.array([S.wrap(th[u, 4][k] - th[u, 1][k]) for u in subs]) for k in ("core", "cue", "locked")}
dcc = np.array([S.wrap((th[u, 4]["core"] - th[u, 4]["cue"]) - (th[u, 1]["core"] - th[u, 1]["cue"])) for u in subs])
cr = []
for k in ("core", "cue", "locked"):
    c, lo, hi = ci(cmean, dphi[k]); cr.append(dict(item=f"timed phase change D4-D1, {k} state (h)", estimate=c, lo=lo, hi=hi,
                                                   resultant=float(S.resultant(dphi[k])), n=len(subs)))
c, lo, hi = ci(cmean, dcc); cr.append(dict(item="timed core-minus-cue phase change (h)", estimate=c, lo=lo, hi=hi,
                                           resultant=float(S.resultant(dcc)), n=len(subs)))
for lab, ph in (("D1 core phase (h clock)", [th[u, 1]["core"] for u in subs]), ("D4 core phase (h clock)", [th[u, 4]["core"] for u in subs])):
    c, lo, hi = ci(cmean, ph); cr.append(dict(item=lab, estimate=c, lo=lo, hi=hi, resultant=float(S.resultant(ph)), n=len(subs)))
for lab, k in (("D1 core cosinor amplitude, mean (state units)", 1), ("D4 core cosinor amplitude, mean (state units)", 4)):
    cr.append(dict(item=lab, estimate=float(np.mean([th[u, k]["core_amp"] for u in subs])), n=len(subs)))


def ccor(a, b):
    a, b = np.asarray(a) / HOURS, np.asarray(b) / HOURS
    sa, sb = np.sin(a - np.angle(np.exp(1j * a).mean())), np.sin(b - np.angle(np.exp(1j * b).mean()))
    return float((sa * sb).sum() / np.sqrt((sa ** 2).sum() * (sb ** 2).sum()))


r = ccor(D["cue"], dcc)
idx = [rb.integers(0, len(subs), len(subs)) for _ in range(NB)]
rbs = np.array([ccor(D["cue"][i], dcc[i]) for i in idx])
rperm = np.array([ccor(D["cue"], dcc[rb.permutation(len(subs))]) for _ in range(NPERM)])
pr = (1 + np.sum(np.abs(rperm) >= abs(r))) / (1 + NPERM)
cr.append(dict(item="circular correlation: untimed d core-cue vs timed core-minus-cue change", estimate=r,
               lo=float(np.nanquantile(rbs, 0.025)), hi=float(np.nanquantile(rbs, 0.975)), p_perm_two_sided=float(pr), n=len(subs)))
cr.append(dict(item="mean wrapped difference untimed minus timed d core-cue (h)", estimate=cmean(S.wrap(D["cue"] - dcc)),
               resultant=float(S.resultant(S.wrap(D["cue"] - dcc))), n=len(subs)))
pd.DataFrame(cr).to_csv(OUT / "descriptive_timed.csv", index=False)
print(tests.to_string()); print(pd.DataFrame(sens).T.to_string()); print(pd.read_csv(OUT / "gene_coverage.csv").to_string())
