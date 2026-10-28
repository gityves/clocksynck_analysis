"""A124 (registered 2026-10-28): does the core-cue offset detect the imposed 12 h misalignment in human muscle (GSE106800)?

Estimator = A115 / A78 (log2 arrays, equal gene weights, one projector per gene set over all primary arrays, states centred
within subject, offset = arg E[c_core c_cue*] per condition via `clocksync.stats.ccorr`). H124.1: |wrapped change of the
core-cue offset| under exact within-subject condition-label swaps (2^11). H124.2: |d core-cue| - |d core-locked|, one-sided,
same swaps. Holm over the two. Seed 124 (bootstrap only; the test is exact).

    python compute_a124.py [--alias-sensitivity]         (run with the sccc python; reads the series matrix and output/probe_map.csv.gz, samples.csv)
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER.parent)); sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A115_A116_psoriasis/scripts"))
import _lib as L  # A115 helpers: states_for (one projector, genes absent dropped), restrict
import clocksync as cs
from clocksync import stats as S
from clocksync.config import limit_threads
from clocksync.genes import canonical

limit_threads(1)
ALIAS_SENS = "--alias-sensitivity" in sys.argv  # post hoc: add the platform-name aliases the project table lacks
EXTRA = {"FAM46B": "TENT5B", "MTERFD3": "MTERF2", "HIST2H2BE": "H2BC21", "C12ORF23": "TMEM263", "SEPT4": "SEPTIN4",
         "C10ORF116": "ADIRF-AS1", "CENPC1": "CENPC", "FAM134B": "RETREG1"}
OUT = HERE.parent / "output"
IN = OUT
if ALIAS_SENS:
    OUT = OUT / "alias_sensitivity"; OUT.mkdir(exist_ok=True)


def canon(s):
    c = canonical(s)
    return EXTRA.get(c, c) if ALIAS_SENS else c
GEO = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/geo/GSE106800")
HOURS = S.HOURS
NB, SEED = 2000, 124
TIMES = ("7 AM (7h)", "7 PM (19h)")  # morning, evening, as recorded

# ---- data ---------------------------------------------------------------------------------------------------------
smp = pd.read_csv(IN / "samples.csv")
smp["cond"] = np.where(smp["experimental condition"].str.contains("misaligned"), "mis", "ctl")
smp["pm"] = (smp["time of sampling"] == TIMES[1]).astype(int)
smp["subject"] = smp.subjectid
pm = pd.read_csv(IN / "probe_map.csv.gz", keep_default_na=False)
pm = pm[pm.category == "main"].copy()
pm["canon"] = pm.symbols_all.map(lambda s: [canon(x) for x in s.split(";") if x])
expr = pd.read_csv(GEO / "GSE106800_series_matrix.txt.gz", sep="\t", comment="!", index_col=0)
expr.index = expr.index.astype(str)
expr = expr[smp.gsm.tolist()]  # columns in the order of samples.csv

sets = {"core": cs.gtex()}
ref = pd.read_csv(PAPER / "results/figure_3/output/gtex_cues/a81b_sets.csv")
for k in ("cue", "locked"):
    x = ref[ref.set == k]
    sets[k] = cs.Waveform(tuple(x.gene), x.amp.to_numpy(), x.acro_h.to_numpy() * 2 * np.pi / 24, None, None, k)


def gene_table():
    """Gene x array table: mean log2 over the gene's main probe sets (those listing it alone if any)."""
    names = sorted({g for w in sets.values() for g in w.genes})
    cols, rows = {}, []
    for g in names:
        c = canon(g)
        m = pm[pm.canon.map(lambda s: c in s)]
        alone = m[m.canon.map(lambda s: s == [c])]
        use = alone if len(alone) else m
        use = use[use.probe.astype(str).isin(expr.index)]
        rows.append(dict(gene=g, canonical=c, n_probesets=len(use), alone=bool(len(alone)), probes=";".join(use.probe.astype(str))))
        if len(use):
            cols[g] = expr.loc[use.probe.astype(str)].mean(axis=0).to_numpy()
    return pd.DataFrame(cols, index=smp.index), pd.DataFrame(rows)


X, cov = gene_table()
for k, w in sets.items():
    cov.loc[cov.gene.isin(w.genes), k] = True
cov.to_csv(OUT / "gene_mapping.csv", index=False)
coverage = {k: (int(cov[cov.gene.isin(w.genes) & (cov.n_probesets > 0)].shape[0]), len(w.genes)) for k, w in sets.items()}

# ---- estimator ----------------------------------------------------------------------------------------------------


def states(rows):
    """Complex states per set for the arrays `rows` (one projector per set over these arrays), centred within subject."""
    x = X.loc[rows]
    sub = smp.loc[rows, "subject"].to_numpy()
    out = {}
    for k, w in sets.items():
        (b, _), n = L.states_for(x.assign(depth=1.0), w, counts=False)
        b = b.copy()
        for s in np.unique(sub):
            b[sub == s] -= b[sub == s].mean(0)
        out[k] = b
        out[k + "_n"] = n
    return out


def off(b1, b2, m):
    z = S.ccorr(b1[m], b2[m])
    return z


def delta(st, cond):
    """(offset mis - ctl) core-cue, core-locked in hours (wrapped), plus the complex values."""
    r = {}
    for k in ("cue", "locked"):
        zc, zm = off(st["core"], st[k], cond == 0), off(st["core"], st[k], cond == 1)
        r[k] = (float(S.wrap((np.angle(zm) - np.angle(zc)) * HOURS)), zc, zm)
    return r


prim_ids = smp.groupby("subject").size()
prim_ids = prim_ids[prim_ids == 4].index.tolist()
rows = smp.index[smp.subject.isin(prim_ids)]
sub = smp.loc[rows, "subject"].to_numpy()
cond0 = (smp.loc[rows, "cond"] == "mis").astype(int).to_numpy()
st = states(rows)
subj_list = sorted(prim_ids)
si = np.array([subj_list.index(s) for s in sub])
n = len(subj_list)

obs = delta(st, cond0)
d_cc, d_cl = obs["cue"][0], obs["locked"][0]
T1, T2 = abs(d_cc), abs(d_cc) - abs(d_cl)

# exact null: all 2^n within-subject swaps (the complement of a swap gives the same |d|, kept for exactness)
null = []
for bits in itertools.product((0, 1), repeat=n):
    c = cond0 ^ np.array(bits)[si]
    r = delta(st, c)
    null.append((r["cue"][0], r["locked"][0]))
null = np.array(null)
a1, a2 = np.abs(null[:, 0]), np.abs(null[:, 0]) - np.abs(null[:, 1])
p1 = float(np.mean(a1 >= T1 - 1e-12)); p2 = float(np.mean(a2 >= T2 - 1e-12))
order = np.argsort([p1, p2]); holm = np.empty(2); run = 0
for rank, i in enumerate(order):
    run = max(run, (2 - rank) * [p1, p2][i]); holm[i] = min(run, 1)
sd_null = float(np.sqrt(-2 * np.log(abs(np.mean(np.exp(1j * null[:, 0] / HOURS))))) * HOURS)  # circular SD of signed null
mde = dict(null_q95_abs=float(np.quantile(a1, 0.95)), null_circ_sd_h=sd_null, mde_80=(1.96 + 0.8416) * sd_null,
           null_q95_h2=float(np.quantile(a2, 0.95)), n_assign=len(null), p_min=1 / len(null))
res = pd.DataFrame([
    dict(test="H124.1", statistic="|d core-cue| (h)", observed=T1, p=p1, p_holm=holm[0], alternative="two-sided", pass_=holm[0] < 0.05),
    dict(test="H124.2", statistic="|d core-cue| - |d core-locked| (h)", observed=T2, p=p2, p_holm=holm[1], alternative="greater", pass_=holm[1] < 0.05)])
res.to_csv(OUT / "tests.csv", index=False)
pd.DataFrame([mde]).to_csv(OUT / "mde.csv", index=False)
pd.DataFrame({"d_core_cue": null[:, 0], "d_core_locked": null[:, 1]}).to_csv(OUT / "null_distribution.csv.gz", index=False)

# ---- descriptives -------------------------------------------------------------------------------------------------
rng = np.random.default_rng(SEED)
boot = []
for _ in range(NB):
    pick = rng.integers(0, n, n)
    idx = np.concatenate([np.where(si == p)[0] for p in pick])
    stb = {k: v[idx] if isinstance(v, np.ndarray) else v for k, v in st.items()}
    r = delta(stb, cond0[idx])
    boot.append((r["cue"][0], r["locked"][0], r["cue"][0] - r["locked"][0]))
boot = np.array(boot)


def ci(v):
    """Percentile interval of wrapped differences around their circular centre."""
    c = S.circmean(v)
    w = S.wrap(np.asarray(v) - c)
    return float(c), float(c + np.quantile(w, 0.025)), float(c + np.quantile(w, 0.975))


def coh(z):
    return float(abs(z))


desc = []
for nm, i, k in (("d core-cue", 0, "cue"), ("d core-locked", 1, "locked")):
    c, lo, hi = ci(boot[:, i])
    desc.append(dict(item=nm, estimate_h=obs[k][0], boot_lo=lo, boot_hi=hi, boot_n=NB))
for k in ("cue", "locked"):
    for c, z in (("ctl", obs[k][1]), ("mis", obs[k][2])):
        desc.append(dict(item=f"offset core-{k} {c}", estimate_h=float(np.angle(z) * HOURS), coherence=coh(z)))
pd.DataFrame(desc).to_csv(OUT / "descriptive_deltas.csv", index=False)

# state difference between the two biopsy times (7 PM minus 7 AM), per condition, per set; raw (uncentred) states
raw = {}
x = X.loc[rows]
for k, w in sets.items():
    (b, _), _ = L.states_for(x.assign(depth=1.0), w, counts=False)
    raw[k] = b[:, 0] + 1j * b[:, 1]
tvec = smp.loc[rows, "pm"].to_numpy()


def tdiff(k, cond, pick=None):
    pick = np.arange(n) if pick is None else pick
    d = []
    for p in pick:
        m = si == p
        e = (cond0 == cond) & m
        d.append(raw[k][e & (tvec == 1)][0] - raw[k][e & (tvec == 0)][0])
    return np.array(d)


rng = np.random.default_rng(SEED + 1)
trows = []
bs = {k: [] for k in ("core", "cue")}
for _ in range(NB):
    pick = rng.integers(0, n, n)
    for k in bs:
        zc, zm = tdiff(k, 0, pick).mean(), tdiff(k, 1, pick).mean()
        bs[k].append(((np.angle(zm) - np.angle(zc)) * HOURS, abs(zc), abs(zm)))
for k in ("core", "cue"):
    zc, zm = tdiff(k, 0).mean(), tdiff(k, 1).mean()
    ph = float(S.wrap((np.angle(zm) - np.angle(zc)) * HOURS))
    c, lo, hi = ci(np.array(bs[k])[:, 0])
    # reference: the GTEx waveform maps clock time t to a state phase; here the state-difference direction is reported relative
    # between conditions (12 h = full reversal of the evening-minus-morning move)
    trows.append(dict(set=k, mean_diff_ctl_abs=float(abs(zc)), mean_diff_mis_abs=float(abs(zm)),
                      dir_ctl_h=float(np.angle(zc) * HOURS), dir_mis_h=float(np.angle(zm) * HOURS),
                      dir_mis_minus_ctl_h=ph, boot_lo=lo, boot_hi=hi,
                      cos_between=float(np.real(zm * np.conj(zc)) / (abs(zm) * abs(zc)))))
    for cnd, z, dd in (("ctl", zc, tdiff(k, 0)), ("mis", zm, tdiff(k, 1))):
        trows[-1][f"t_amp_over_se_{cnd}"] = float(abs(z) / (np.sqrt(np.var(dd.real, ddof=1) + np.var(dd.imag, ddof=1)) / np.sqrt(n)))
pd.DataFrame(trows).to_csv(OUT / "descriptive_time_state_difference.csv", index=False)

# sensitivity (descriptive): all 12 subjects, projector over all 47 arrays, no test
rows_all = smp.index
sta = states(rows_all)
ca = (smp.loc[rows_all, "cond"] == "mis").astype(int).to_numpy()
ra = delta(sta, ca)
pd.DataFrame([dict(n_arrays=len(rows_all), n_subjects=smp.subject.nunique(), d_core_cue_h=ra["cue"][0], d_core_locked_h=ra["locked"][0])]
             ).to_csv(OUT / "sensitivity_12_subjects.csv", index=False)

pd.DataFrame([dict(set=k, genes_in_set=coverage[k][1], genes_mapped=coverage[k][0], used_in_projector=int(st[k + "_n"])) for k in sets]
             ).to_csv(OUT / "gene_coverage.csv", index=False)
print(res.to_string()); print(mde); print(pd.read_csv(OUT / "gene_coverage.csv").to_string())
