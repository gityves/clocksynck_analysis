"""A255: do two different tissues of a living person share a clock state, beyond what the person's genome gives?

Registered in paper/docs/registered_experiments.md (A255). MuTHER twins (fat, skin, lymphoblastoid cell line) and three
surgical two-tissue cohorts. State: clock-plane state (GTEx waveform, unweighted projector), centred within tissue;
pairing r = |sum za conj(zb)| / sqrt(sum |za|^2 sum |zb|^2) over the persons with both; null: persons permuted within tissue.

Stages
  power  persons per pair, gates, null SD and MDE of each r (permuted persons only) -> output/{power,gates}.csv
  score  (after the power note is dated) everything registered -> output/score_{tests,twins,descriptive}.csv
  `--placebo`: persons shuffled once within each tissue before anything else (code test; output/placebo/).

    python compute_a255.py --stage power|score [--smoke] [--placebo]

Choices not fixed by the registration
- Input: `A254_paired_lesions/output/prep/` (`prep_lesions_tissues.py`): MuTHER matrices per tissue (rows = women) and
  `muther_persons.csv`; surgical cohorts as expression + sheet. Arrays: deposited log2 values x ln 2; variance-stabilised counts
  (GSE133979) and two-colour log ratios (GSE24883) are used as deposited (only centred states enter).
- A clock gene absent from a cohort is dropped from the waveform there.
- Companion sets: drawn from the mean over the two tissues (MuTHER: fat and skin) of the mean log expression; the same gene
  sets in all tissues of a cohort; set j is pooled across the surgical cohorts.
- H255.2: women resampled with replacement (2,000); each r on the resampled women having both tissues.
- H255.4: one twin pair = two women listed as each other's co-twin; the order within a pair is randomised once (seed);
  null: the second twins permuted over pairs (10,000); the difference to the cell line is bootstrapped over twin pairs.
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
PREP = PAPER / "exploratory/A254_paired_lesions/output/prep"
SEED, N_PERM, N_SETS, N_BOOT, MDE_K = 255, 10_000, 200, 2000, 2.49
WF = cs.gtex()
HOURS = 24 / (2 * np.pi)
SURG = ["GSE73439", "GSE133979", "GSE24883"]


def say(*a):
    import time
    print(time.strftime("%F %T"), *a, flush=True)


def zstates(y, wf):
    g = list(wf.genes)
    L = P.projector(y[g], wf, counts=False)
    b = y[g].to_numpy(float) @ L.T
    return b[:, 0] + 1j * b[:, 1]


class Cohort:
    """tissues: name -> DataFrame (persons x genes, log scale); all aligned to `ids`."""

    def __init__(self, name, tissues, placebo=False):
        self.name, self.t = name, tissues
        self.ids = sorted(set().union(*[set(v.index) for v in tissues.values()]))
        common = set.intersection(*[set(v.columns) for v in tissues.values()])
        self.wf = WF.subset([g for g in WF.genes if g in common])
        if placebo:
            rng = np.random.default_rng([SEED, 9, sum(map(ord, name))])
            for k, v in tissues.items():
                tissues[k] = pd.DataFrame(v.to_numpy()[rng.permutation(len(v))], index=v.index, columns=v.columns)
        two = list(tissues)[:2]
        lm = sum(tissues[k].mean() for k in two) / 2
        self.lm = lm[[c in common for c in lm.index]].dropna()
        self.lm = self.lm[~self.lm.index.duplicated()]

    def Z(self, wf=None):
        wf = wf or self.wf
        out = {}
        for k, v in self.t.items():
            z = pd.Series(zstates(v, wf), index=v.index).reindex(self.ids).to_numpy()
            out[k] = z - np.nanmean(z)
        return out

    def pseudo(self, kind, n, tag):
        return P.pseudo_waveforms(self.lm, n, np.random.default_rng([SEED, tag, sum(map(ord, self.name))]), self.wf, kind)

    def gate(self, k, i):
        return P.reference_gate(self.t[k][list(self.wf.genes)].to_numpy(float), self.wf, 200, np.random.default_rng([SEED, 1, i]))


def r_pair(za, zb):
    ok = ~(np.isnan(za) | np.isnan(zb))
    a, b = za[ok], zb[ok]
    c = np.sum(a * np.conj(b))
    return float(abs(c) / np.sqrt(np.sum(abs(a) ** 2) * np.sum(abs(b) ** 2))), c, int(ok.sum())


def perm_null(za, zb, n_perm, rng):
    ok = ~(np.isnan(za) | np.isnan(zb))
    a, b = za[ok], zb[ok]
    den = np.sqrt(np.sum(abs(a) ** 2) * np.sum(abs(b) ** 2))
    perms = np.argsort(rng.random((n_perm, len(a))), axis=1)
    return np.abs((a[None, :] * np.conj(b[perms])).sum(1)) / den


def p_up(nul, obs):
    return float((1 + np.sum(np.asarray(nul) >= obs)) / (1 + len(nul)))


def load_muther(placebo):
    t = {k: pd.read_parquet(PREP / f"muther_{k}_expr.parquet").astype(float) * np.log(2.0) for k in ("fat", "skin", "lcl")}
    return Cohort("muther", t, placebo), pd.read_csv(PREP / "muther_persons.csv", dtype=str)


def load_surg(acc, placebo):
    x = pd.read_parquet(PREP / f"{acc}_expr.parquet").astype(float)
    s = pd.read_csv(PREP / f"{acc}_sheet.csv", dtype={"person": str}).set_index("sample")
    x = x.loc[s.index]
    if (s.kind == "counts").all():
        x = np.log(x + 1.0)
    elif not s.scale_note.astype(str).str.contains("ratio|stabil", case=False).any():
        x = x * np.log(2.0)
    t = {}
    for tis in sorted(s.tissue.unique()):
        v = x[(s.tissue == tis).to_numpy()]
        v.index = s.person[v.index].to_numpy()
        t[tis] = v[~v.index.duplicated()].dropna(axis=1)
    return Cohort(acc, t, placebo)


def companions(c, pairs, n_sets):
    out = {}
    for kind, tag in (("matched", 3), ("responsive", 4)):
        vals = {p: [] for p in pairs}
        for w in c.pseudo(kind, n_sets, tag):
            Zp = c.Z(w)
            for a, b in pairs:
                vals[(a, b)].append(r_pair(Zp[a], Zp[b])[0])
        out[kind] = {p: np.array(v) for p, v in vals.items()}
    return out


def run(args):
    n_perm, n_sets, n_boot = (300, 5, 100) if args.smoke else (N_PERM, N_SETS, N_BOOT)
    M, persons = load_muther(args.placebo)
    S = {a: load_surg(a, args.placebo) for a in SURG}
    Z = M.Z()
    mp = [("fat", "skin"), ("fat", "lcl"), ("skin", "lcl")]
    gates = [dict(cohort="muther", tissue=k, persons=len(M.t[k]), clock_genes=len(M.wf.genes), **M.gate(k, i)) for i, k in enumerate(M.t)]
    for a, c in S.items():
        gates += [dict(cohort=a, tissue=k, persons=len(c.t[k]), clock_genes=len(c.wf.genes), **c.gate(k, 10 + i)) for i, k in enumerate(c.t)]
    G = pd.DataFrame(gates)
    G.to_csv(OUT / "gates.csv", index=False)
    rng = np.random.default_rng([SEED, 2])
    nul = {p: perm_null(Z[p[0]], Z[p[1]], n_perm, rng) for p in mp}
    ZS = {a: c.Z() for a, c in S.items()}
    sp = {a: tuple(c.t)[:2] for a, c in S.items()}
    nulS = {a: perm_null(ZS[a][sp[a][0]], ZS[a][sp[a][1]], n_perm, rng) for a in S}
    wS = np.array([r_pair(ZS[a][sp[a][0]], ZS[a][sp[a][1]])[2] for a in S], float)
    pooledS = sum(w * nulS[a] for w, a in zip(wS / wS.sum(), S))
    power = [dict(test=f"muther {a}-{b}", persons=r_pair(Z[a], Z[b])[2], null_mean=float(nul[(a, b)].mean()), null_sd=float(nul[(a, b)].std()),
                  mde=MDE_K * float(nul[(a, b)].std())) for a, b in mp]
    power += [dict(test=f"{a} {sp[a][0]}-{sp[a][1]}", persons=int(w), null_mean=float(nulS[a].mean()), null_sd=float(nulS[a].std()), mde=MDE_K * float(nulS[a].std()))
              for a, w in zip(S, wS)]
    power.append(dict(test="surgical cohorts pooled", persons=int(wS.sum()), null_mean=float(pooledS.mean()), null_sd=float(pooledS.std()), mde=MDE_K * float(pooledS.std())))
    Pw = pd.DataFrame(power)
    Pw.to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 250)
    print(G.round(4).to_string(), Pw.round(4).to_string(), sep="\n", flush=True)
    if args.stage == "power":
        return

    comp = companions(M, mp, n_sets)
    tests = []
    for a, b in mp:
        r, c, n = r_pair(Z[a], Z[b])
        row = dict(test=f"muther {a}-{b}", persons=n, r=r, null_mean=float(nul[(a, b)].mean()), null_sd=float(nul[(a, b)].std()), p=p_up(nul[(a, b)], r),
                   offset_h=float(np.angle(c) * HOURS))
        for kind in comp:
            v = comp[kind][(a, b)]
            row.update({f"{kind}_median": float(np.median(v)), f"{kind}_q95": float(np.quantile(v, .95)), f"p_{kind}": p_up(v, r)})
        row["clock_specific"] = bool(r > row["matched_q95"] and r > row["responsive_q95"])
        tests.append(row)
    compS = {a: companions(S[a], [sp[a]], n_sets) for a in S}
    obsS = []
    for a in S:
        r, c, n = r_pair(ZS[a][sp[a][0]], ZS[a][sp[a][1]])
        obsS.append(r)
        row = dict(test=f"{a} {sp[a][0]}-{sp[a][1]}", persons=n, r=r, null_mean=float(nulS[a].mean()), null_sd=float(nulS[a].std()), p=p_up(nulS[a], r),
                   offset_h=float(np.angle(c) * HOURS))
        for kind in compS[a]:
            v = compS[a][kind][sp[a]]
            row.update({f"{kind}_median": float(np.median(v)), f"{kind}_q95": float(np.quantile(v, .95)), f"p_{kind}": p_up(v, r)})
        row["clock_specific"] = bool(r > row["matched_q95"] and r > row["responsive_q95"])
        tests.append(row)
    ro = float(np.sum(wS / wS.sum() * obsS))
    row = dict(test="surgical cohorts pooled", persons=int(wS.sum()), r=ro, null_mean=float(pooledS.mean()), null_sd=float(pooledS.std()), p=p_up(pooledS, ro))
    for kind in ("matched", "responsive"):
        v = sum(w * compS[a][kind][sp[a]] for w, a in zip(wS / wS.sum(), S))
        row.update({f"{kind}_median": float(np.median(v)), f"{kind}_q95": float(np.quantile(v, .95)), f"p_{kind}": p_up(v, ro)})
    row["clock_specific"] = bool(ro > row["matched_q95"] and ro > row["responsive_q95"])
    tests.append(row)
    T = pd.DataFrame(tests)

    desc = []
    rng = np.random.default_rng([SEED, 5])
    n = len(M.ids)
    diff = lambda Zx: r_pair(Zx["fat"], Zx["skin"])[0] - 0.5 * (r_pair(Zx["fat"], Zx["lcl"])[0] + r_pair(Zx["skin"], Zx["lcl"])[0])  # noqa: E731
    bs = np.array([diff({k: v[i] for k, v in Z.items()}) for i in (rng.integers(0, n, n) for _ in range(n_boot))])
    desc.append(dict(item="H255.2 r(fat, skin) - mean r(tissue, cell line)", value=diff(Z), lo90=float(np.quantile(bs, .05)), hi90=float(np.quantile(bs, .95))))
    Z2 = M.Z(M.wf.drop(["PER1", "NR1D1"]))
    desc += [dict(item=f"r without PER1 and NR1D1, muther {a}-{b}", value=r_pair(Z2[a], Z2[b])[0]) for a, b in mp]
    c = r_pair(Z["fat"], Z["skin"])[1]
    ok = np.flatnonzero(~(np.isnan(Z["fat"]) | np.isnan(Z["skin"])))
    off = float(np.angle(c) * HOURS)
    bo = np.array([np.angle(np.sum(Z["fat"][i] * np.conj(Z["skin"][i]))) * HOURS for i in (rng.choice(ok, len(ok)) for _ in range(n_boot))])
    bo = (bo - off + 12) % 24 - 12 + off
    desc.append(dict(item="fat minus skin offset (h)", value=off, lo90=float(np.quantile(bo, .05)), hi90=float(np.quantile(bo, .95))))

    # co-twins
    pid = pd.Index(M.ids)
    pr = persons[persons.cotwin.isin(pid) & persons.person.isin(pid)]
    pairs = sorted({tuple(sorted((a, b))) for a, b in zip(pr.person, pr.cotwin) if a != b})
    flip = np.random.default_rng([SEED, 6]).random(len(pairs)) < 0.5
    i1 = pid.get_indexer([b if f else a for (a, b), f in zip(pairs, flip)])
    i2 = pid.get_indexer([a if f else b for (a, b), f in zip(pairs, flip)])
    zyg = persons.set_index("person").zygosity.reindex([a for a, _ in pairs]).to_numpy()
    tw = []
    rng = np.random.default_rng([SEED, 7])
    rt = lambda k, sel: r_pair(Z[k][i1[sel]], Z[k][i2[sel]])  # noqa: E731
    allp = np.arange(len(pairs))
    for k in ("fat", "skin", "lcl"):
        r, _, npair = rt(k, allp)
        nl = perm_null(Z[k][i1], Z[k][i2], n_perm, rng)
        row = dict(tissue=k, twin_pairs=npair, r=r, null_mean=float(nl.mean()), null_sd=float(nl.std()), p=p_up(nl, r))
        for zy in ("MZ", "DZ"):
            sel = np.flatnonzero(zyg == zy)
            if len(sel) >= 20:
                row[f"r_{zy}"], _, row[f"pairs_{zy}"] = rt(k, sel)
        if k != "lcl":
            d = r - rt("lcl", allp)[0]
            bsd = np.array([rt(k, s)[0] - rt("lcl", s)[0] for s in (rng.integers(0, len(pairs), len(pairs)) for _ in range(n_boot))])
            row.update(minus_cell_line=d, lo90=float(np.quantile(bsd, .05)), hi90=float(np.quantile(bsd, .95)))
        tw.append(row)
    TW = pd.DataFrame(tw)
    T.to_csv(OUT / "score_tests.csv", index=False)
    TW.to_csv(OUT / "score_twins.csv", index=False)
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    print(T.round(4).T.to_string(), TW.round(4).T.to_string(), pd.DataFrame(desc).round(4).to_string(), sep="\n", flush=True)


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
