"""A251: do skin compartments of distant body sites of a living person share one clock state? (Wu 2020, GSE139301 + GSE139300)

Registered in paper/docs/registered_experiments.md (A251). Six compartments per woman (epidermis / dermis x forearm /
buttock / cheek), each hybridised in its own batch; untimed. State: clock-plane state of the 14 genes (GTEx waveform,
unweighted projector), centred within compartment. Pair statistic r = |sum za conj(zb)| / sqrt(sum |za|^2 sum |zb|^2);
families pooled with weights = women; null: women permuted independently within each compartment.

Stages
  power  design, gate per compartment, null SD and MDE of each family (permuted women only) -> output/{power,gates}.csv
  score  (after the power note is dated) everything registered -> output/score_{tests,pairs,offsets,timed,descriptive}.csv
  `--placebo`: women shuffled once within each compartment before anything else (code test; output/placebo/).

    python compute_a251.py --stage power|score [--smoke] [--placebo]

Choices not fixed by the registration
- Matrix values are used as log2 intensities (GEO: RMA); only centred states enter the statistics, so the base of the
  logarithm changes no pairing, offset or p.
- Gene symbols: column 15 of GPL13667 (first symbol of a `///` list), canonical names (`clocksync.genes.canonical`).
- Pseudo-clock sets: drawn once from the mean over the two layers of the women's mean log intensity per gene, the same
  gene sets in all compartments; a set is dropped if a gene is missing in a layer.
- Family p: share of permutations with pooled r >= observed. Companion p: share of sets >= observed.
- H251.3: percentile bootstrap of the women of each pair (2,000), offsets wrapped to +/- 12 h around the estimate.
- H251.4: men resampled independently per layer (the dermis and epidermis numbers of a man are not linked in GEO).
- P251.0: a = mean of z exp(-i w t) over samples, z centred within man; null: the four times permuted within man.
"""

from __future__ import annotations

import gzip
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parent / "output"
SRC = DATA / "skin_wu2020_gse139301"
SHEET = PAPER / "exploratory/scouts_2026-11-12/output/wu2020_GSE139300_GSE139301_samples.csv"
SEED, N_PERM, N_SETS, N_BOOT, MDE_K = 251, 10_000, 200, 2000, 2.49
WF = cs.gtex()
LAYERS, SITES = ("epidermis", "dermis"), ("forearm", "buttock", "cheek")
COMPS = [f"{l}_{s}" for l in LAYERS for s in SITES]
W = 2 * np.pi / 24
HOURS = 1 / W


def families():
    same_layer = [(f"{l}_{a}", f"{l}_{b}") for l in LAYERS for a, b in combinations(SITES, 2)]
    same_site = [(f"epidermis_{s}", f"dermis_{s}") for s in SITES]
    cross = [(f"epidermis_{a}", f"dermis_{b}") for a in SITES for b in SITES if a != b]
    return {"H251.1 same layer, different sites": same_layer, "H251.2a layers, same site": same_site,
            "H251.2b layers, different sites": cross}


def load():
    """layer -> (samples x genes log2 matrix, sheet rows), genes collapsed by the registered rule."""
    sheet = pd.read_csv(SHEET).set_index("sample_id")
    sym = {}
    with open(SRC / "GPL13667_full.txt", errors="replace") as f:
        for line in f:
            if line[0] in "#^!":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > 14 and p[14] and p[14] != "---" and p[0] != "ID":
                sym[p[0]] = canonical(p[14].split("///")[0].strip())
    out = {}
    for layer, gse in (("epidermis", "GSE139301"), ("dermis", "GSE139300")):
        with gzip.open(SRC / f"{gse}_series_matrix.txt.gz", "rt") as f:
            skip = next(i for i, l in enumerate(f) if l.startswith("!series_matrix_table_begin")) + 1
        x = pd.read_csv(SRC / f"{gse}_series_matrix.txt.gz", sep="\t", skiprows=skip, index_col=0, comment="!")
        x = x.loc[x.index.isin(sym)].astype(np.float32)
        s = sheet.loc[x.columns]
        women = s.index[s.arm != "timed_men"]
        g = pd.Series({i: sym[i] for i in x.index})
        mean = x[women].mean(1)
        specific = ~g.index.str.contains("_s_|_x_")
        has_spec = pd.Series(specific, index=g.index).groupby(g).transform("any")
        cand = mean[(specific | ~has_spec.to_numpy())]
        best = cand.groupby(g[cand.index]).idxmax()
        m = x.loc[best.to_numpy()].T
        m.columns = best.index
        out[layer] = (m, s)
    return out


def z_of(m, wf):
    g = list(wf.genes)
    L = P.projector(m[g], wf, counts=False)
    b, _ = P.states(m[g], L, wf, depth=None)
    return b[:, 0] + 1j * b[:, 1]


def women_states(data, wf, ids):
    """compartment -> complex states over `ids` (NaN where missing), centred within compartment."""
    Z = {}
    for layer in LAYERS:
        m, s = data[layer]
        if any(g not in m.columns for g in wf.genes):
            return None
        z = pd.Series(z_of(m, wf), index=m.index)
        for site in SITES:
            k = s.index[(s.arm != "timed_men") & (s.site == site)]
            v = pd.Series(z[k].to_numpy(), index=s.loc[k, "subject_id"]).reindex(ids).to_numpy()
            Z[f"{layer}_{site}"] = v - np.nanmean(v)
    return Z


def r_pair(za, zb):
    ok = ~(np.isnan(za) | np.isnan(zb))
    a, b = za[ok], zb[ok]
    c = np.sum(a * np.conj(b))
    return abs(c) / np.sqrt(np.sum(abs(a) ** 2) * np.sum(abs(b) ** 2)), c, int(ok.sum())


def pooled(Z, pairs):
    v = [r_pair(Z[a], Z[b]) for a, b in pairs]
    w = np.array([x[2] for x in v], float)
    return float(np.sum(w * [x[0] for x in v]) / w.sum())


def null(Z, fams, n_perm, rng):
    n = len(next(iter(Z.values())))
    out = {k: np.empty(n_perm) for k in fams}
    for i in range(n_perm):
        Zp = {c: v[rng.permutation(n)] for c, v in Z.items()}
        for k, pairs in fams.items():
            out[k][i] = pooled(Zp, pairs)
    return out


def p_up(nul, obs):
    return float((1 + np.sum(np.asarray(nul) >= obs)) / (1 + len(nul)))


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def offsets(Z, pairs, n_boot, rng):
    rows = []
    for a, b in pairs:
        ok = np.flatnonzero(~(np.isnan(Z[a]) | np.isnan(Z[b])))
        za, zb = Z[a][ok], Z[b][ok]
        est = np.angle(np.sum(za * np.conj(zb))) * HOURS
        bs = np.array([np.angle(np.sum(za[i] * np.conj(zb[i]))) * HOURS for i in (rng.integers(0, len(ok), len(ok)) for _ in range(n_boot))])
        d = wrap(bs - est) + est
        rows.append(dict(a=a, b=b, women=len(ok), offset_h=float(est), lo90=float(np.quantile(d, .05)), hi90=float(np.quantile(d, .95))))
    return pd.DataFrame(rows)


def timed(data, wf, n_perm, n_boot, rng):
    """Per layer: complex cosinor of the men's states on clock time; null: times permuted within man."""
    res, a_boot = [], {}
    for layer in LAYERS:
        m, s = data[layer]
        k = s.index[s.arm == "timed_men"]
        z = pd.Series(z_of(m, wf), index=m.index)[k].to_numpy()
        man = s.loc[k, "subject_id"].to_numpy()
        t = s.loc[k, "collection_time_hhmm"].to_numpy(float) / 100.0
        z = z - pd.Series(z).groupby(man).transform("mean").to_numpy()
        a = np.mean(z * np.exp(-1j * W * t))
        idx = [np.flatnonzero(man == u) for u in np.unique(man)]
        nul = np.empty(n_perm)
        for i in range(n_perm):
            tp = t.copy()
            for j in idx:
                tp[j] = t[rng.permutation(j)]
            nul[i] = abs(np.mean(z * np.exp(-1j * W * tp)))
        a_boot[layer] = np.array([np.mean(np.concatenate([z[idx[u]] * np.exp(-1j * W * t[idx[u]]) for u in rng.integers(0, len(idx), len(idx))]))
                                  for _ in range(n_boot)])
        res.append(dict(layer=layer, men=len(idx), samples=len(k), modulus=float(abs(a)), null_mean=float(nul.mean()), null_q95=float(np.quantile(nul, .95)),
                        p=p_up(nul, abs(a)), state_minus_clock_h=float(np.angle(a) * HOURS), _a=a))
    T = pd.DataFrame(res)
    a_e, a_d = T.loc[T.layer == "epidermis", "_a"].iloc[0], T.loc[T.layer == "dermis", "_a"].iloc[0]
    est = float(np.angle(a_d * np.conj(a_e)) * HOURS)
    d = wrap(np.angle(a_boot["dermis"] * np.conj(a_boot["epidermis"])) * HOURS - est) + est
    return T.drop(columns="_a"), dict(men_offset_h=est, men_lo90=float(np.quantile(d, .05)), men_hi90=float(np.quantile(d, .95)))


def run(args):
    data = load()
    n_perm, n_sets, n_boot = (200, 5, 100) if args.smoke else (N_PERM, N_SETS, N_BOOT)
    ids = sorted(set(data["epidermis"][1].query("arm != 'timed_men'").subject_id) | set(data["dermis"][1].query("arm != 'timed_men'").subject_id))
    rng = np.random.default_rng([SEED, 0])
    if args.placebo:
        for layer in LAYERS:
            m, s = data[layer]
            v = m.to_numpy().copy()
            for site in SITES:
                k = np.flatnonzero(((s.arm != "timed_men") & (s.site == site)).to_numpy())
                v[k] = v[rng.permutation(k)]
            data[layer] = (pd.DataFrame(v, index=m.index, columns=m.columns), s)
    fams = families()
    Z = women_states(data, WF, ids)
    gates = []
    for layer in LAYERS:
        m, s = data[layer]
        for site in SITES:
            k = s.index[(s.arm != "timed_men") & (s.site == site)]
            g = P.reference_gate(m.loc[k, list(WF.genes)].to_numpy(float) * np.log(2), WF, 200, np.random.default_rng([SEED, 1, len(gates)]))
            gates.append(dict(compartment=f"{layer}_{site}", women=len(k), **g))
    G = pd.DataFrame(gates)
    G.to_csv(OUT / "gates.csv", index=False)
    nul = null(Z, fams, n_perm, np.random.default_rng([SEED, 2]))
    power = pd.DataFrame([dict(family=k, pairs=len(p), women_min=min(r_pair(Z[a], Z[b])[2] for a, b in p), null_mean=float(nul[k].mean()),
                               null_sd=float(nul[k].std()), null_q95=float(np.quantile(nul[k], .95)), mde=MDE_K * float(nul[k].std())) for k, p in fams.items()])
    power["gated_compartments"] = int((G.p < 0.05).sum())
    power["scored"] = bool((G.p < 0.05).sum() >= 4 and power.mde.iloc[0] <= 0.10)
    power.to_csv(OUT / "power.csv", index=False)
    print(G.round(4).to_string(), power.round(4).to_string(), sep="\n", flush=True)
    if args.stage == "power":
        return

    lm = (sum(data[l][0].loc[data[l][1].arm != "timed_men"].mean() for l in LAYERS) / 2 * np.log(2)).dropna()
    lm = lm[~lm.index.duplicated()]
    comp = {}
    for kind, tag in (("matched", 3), ("responsive", 4)):
        vals = {k: [] for k in fams}
        for w in P.pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, tag]), WF, kind):
            Zp = women_states(data, w, ids)
            if Zp is None:
                continue
            for k, p in fams.items():
                vals[k].append(pooled(Zp, p))
        comp[kind] = {k: np.array(v) for k, v in vals.items()}
    tests = []
    for k, p in fams.items():
        obs = pooled(Z, p)
        row = dict(family=k, pairs=len(p), r=obs, null_mean=float(nul[k].mean()), null_sd=float(nul[k].std()), p=p_up(nul[k], obs), mde=MDE_K * float(nul[k].std()))
        for kind in comp:
            v = comp[kind][k]
            row.update({f"{kind}_median": float(np.median(v)), f"{kind}_q95": float(np.quantile(v, .95)), f"p_{kind}": p_up(v, obs), f"{kind}_sets": len(v)})
        row["clock_specific"] = bool(obs > row["matched_q95"] and obs > row["responsive_q95"])
        tests.append(row)
    T = pd.DataFrame(tests)
    allpairs = list(combinations(COMPS, 2))
    wf2 = WF.drop(["PER1", "NR1D1"])
    Z2 = women_states(data, wf2, ids)
    prow = []
    for a, b in allpairs:
        r, c, n = r_pair(Z[a], Z[b])
        prow.append(dict(a=a, b=b, women=n, r=float(r), r_without_PER1_NR1D1=float(r_pair(Z2[a], Z2[b])[0]),
                         matched_median=float(np.median([0]))))
    Pr = pd.DataFrame(prow).drop(columns="matched_median")
    O = offsets(Z, allpairs, n_boot, np.random.default_rng([SEED, 5]))
    O["within_1.5h"] = (O.lo90 > -1.5) & (O.hi90 < 1.5)
    Tm, men = timed(data, WF, n_perm, n_boot, np.random.default_rng([SEED, 6]))
    desc = [dict(item=f"clock-plane variance, {c}", value=float(np.nanmean(abs(Z[c]) ** 2))) for c in COMPS]
    desc += [dict(item=f"pooled r without PER1 and NR1D1, {k}", value=pooled(Z2, p)) for k, p in fams.items()]
    fo = O[(O.a == "epidermis_forearm") & (O.b == "dermis_forearm")].iloc[0]
    desc += [dict(item="women forearm dermis minus epidermis offset (h)", value=-fo.offset_h, lo90=-fo.hi90, hi90=-fo.lo90),
             dict(item="timed men forearm dermis minus epidermis offset (h)", value=men["men_offset_h"], lo90=men["men_lo90"], hi90=men["men_hi90"])]
    T.to_csv(OUT / "score_tests.csv", index=False)
    Pr.to_csv(OUT / "score_pairs.csv", index=False)
    O.to_csv(OUT / "score_offsets.csv", index=False)
    Tm.to_csv(OUT / "score_timed.csv", index=False)
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    pd.set_option("display.width", 250)
    print(T.round(4).T.to_string(), Pr.round(3).to_string(), O.round(2).to_string(), Tm.round(4).to_string(), pd.DataFrame(desc).round(3).to_string(), sep="\n", flush=True)


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
