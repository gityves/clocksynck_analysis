"""A337 scoring (registered in `paper/docs/registered_experiments.md`): does the A328 lead statistic, blind to time,
recover the pre-mRNA lead of tissues whose clock runs (baboon, Greenwell mouse liver; Perrin muscle when delivered)?

Per compartment, spliced (`counts_exon`) and unspliced (`counts_intron`; `counts_intron_direct` as sensitivity) log
fractions are residualised gene by gene on the covariates (none for baboon; regimen for Greenwell; person for Perrin),
centred and projected with one projector built from the spliced counts (`clocksync.paired.projector`, GTEx waveform;
states turn counter-clockwise with clock time, so arg z is a phase in clock time).
- Untimed (A328): c = sum z_u conj(z_s); gate kappa = |c| / sqrt(sum|z_u|^2 sum|z_s|^2) against permutations of the
  unspliced rows within groups; lead delta_U = arg c (h, > 0: unspliced ahead); 90 % interval from bootstrap
  resamples of samples within groups; 200 matched and 200 responsive companions (own projector, 1,000 permutations).
- Timed: beta = sum z exp(-i 2 pi t / 24); delta_T = arg(beta_u conj(beta_s)); rhythm check |beta_s| against
  permutations of t within groups.
- H337.1: rhythmic, gate (Holm over rhythmic compartments) and delta_U interval within (0, 6 h). H337.2: TOST of
  delta_U - delta_T within +-1.5 h (paired bootstrap; Holm). H337.3: H337.1 and delta_U above the responsive q95 of
  the companions whose own gate passes.
`--check` replaces the unspliced counts with thinned spliced counts and shuffles the times (no real statistic is
computed; expected delta near 0) to test the code.

    python score_a337.py [--sets baboon greenwell | perrin] [--jobs 8] [--check]
Outputs (`../output/<sfx>/`): compartments.csv, verdicts.csv, companions.csv, genes.csv.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import projector, pseudo_waveforms  # noqa: E402
from common import development as dev  # noqa: E402

DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
SRC = DATA / "sccc_mapping"
OUT = HERE.parents[1] / "output"
SEED, N_PERM, N_PERM_PSEUDO, N_BOOT, N_PSEUDO = 337, 10000, 1000, 2000, 200
W = 2 * np.pi / 24
WF = cs.gtex()
WAVEFORMS = {"all": WF, "core": WF.drop(["PER1", "PER2", "NR1D1"])}
GREENWELL_EXCLUDED = ("SRR7747184", "SRR7747208")  # exon counts match each other's processed column (A337)
REGIMEN = {"ad-libitum": "AL", "night-restricted": "NR", "did not adapt": "DR", "arrhythmic": "AR"}
S = {}


def by_symbol(path, cols, symbols) -> pd.DataFrame:
    t = pd.read_csv(path, sep="\t", index_col=0)
    t = t[cols]
    return t.groupby(symbols).sum().T


def load(name) -> list[dict]:
    """Compartments of one set: {name, primary, spliced, unspliced, direct (samples x symbols), obs (t, group)}."""
    d = SRC / {"baboon": "mure2018_baboon", "greenwell": "greenwell2019_mouse_liver",
               "perrin": "perrin2018_muscle"}[name]
    s = pd.read_csv(d / "samples.tsv", sep="\t")
    ann = pd.read_csv(d / "counts_exon.tsv.gz", sep="\t", index_col=0, usecols=[0, 1]).symbol
    stable = pd.Series(ann.index.str.split(".").str[0], index=ann.index)
    if name == "baboon":  # missing symbols from the Ensembl 90 table of the authors, by stable id
        e90 = pd.read_csv(DATA / "baboon/GSE98965/GSE98965_baboon_tissue_expression_FPKM.csv.gz",
                          usecols=["EnsemblID", "Symbol"]).set_index("EnsemblID").Symbol
        ann = ann.fillna(stable.map(e90))
    sym = ann.fillna(stable).astype(str).map(canonical).to_numpy()
    if name == "baboon":
        s["comp"], s["t"], s["group"] = s.char_tissue, s.char_time.str[2:].astype(float), "all"
        units = [(c, c != "Adrenal medulla", s.comp == c) for c in sorted(s.comp.unique())]  # medulla descriptive
    elif name == "greenwell":
        s = s[~s.run_accession.isin(GREENWELL_EXCLUDED)].copy()
        src = s.geo_source_name_ch1.str.lower()
        s["group"] = [next(v for k, v in REGIMEN.items() if k in x) for x in src]
        s["t"] = s.geo_title.str[-2:].astype(float)
        units = [("liver AL+NR", True, s.group.isin(["AL", "NR"])), ("liver AR", False, s.group == "AR"),
                 ("liver DR", False, s.group == "DR")]
    else:  # perrin: columns to be fixed from samples.tsv at delivery, before scoring (person, clock time)
        raise NotImplementedError("Perrin: map person and clock time from samples.tsv after delivery")
    runs = s.run_accession.tolist()
    tabs = {k: by_symbol(d / f"counts_{f}.tsv.gz", runs, sym)
            for k, f in (("spliced", "exon"), ("unspliced", "intron"), ("direct", "intron_direct"))}
    out = []
    for comp, primary, m in units:
        r = s.run_accession[m.to_numpy()].tolist()
        out.append(dict(name=f"{name}:{comp}", primary=primary, obs=s.set_index("run_accession").loc[r, ["t", "group"]],
                        **{k: v.loc[r] for k, v in tabs.items()}))
    return out


def centred_states(counts, obs, wf, L) -> np.ndarray:
    """Complex states of the rows: log fractions residualised on the group indicator, centred, projected."""
    g = list(wf.genes)
    y = np.log((counts[g].to_numpy(float) + 0.5) / counts.sum(1).to_numpy(float)[:, None])
    Z = dev.design(obs, age=None, factors=("group",))
    y = dev.residual(y, Z)
    b = (y - y.mean(0)) @ L.T
    return b[:, 0] + 1j * b[:, 1]


def wrap(h):
    return (np.asarray(h) + 12) % 24 - 12


def lead(zu, zs) -> float:
    return float(np.angle(np.sum(zu * np.conj(zs))) * 24 / (2 * np.pi))


def kappa(zu, zs) -> float:
    return float(np.abs(np.sum(zu * np.conj(zs))) / np.sqrt(np.sum(np.abs(zu) ** 2) * np.sum(np.abs(zs) ** 2)))


def within(groups, rng, n):
    """n permutations (rows) of the indices, within groups."""
    idx = np.tile(np.arange(len(groups)), (n, 1))
    for g in np.unique(groups):
        i = np.flatnonzero(groups == g)
        idx[:, i] = i[rng.permuted(np.tile(np.arange(len(i)), (n, 1)), axis=1)]
    return idx


def boot(groups, rng, n):
    idx = np.empty((n, len(groups)), int)
    for g in np.unique(groups):
        i = np.flatnonzero(groups == g)
        idx[:, i] = rng.choice(i, (n, len(i)))
    return idx


def gate_p(zu, zs, groups, rng, n) -> float:
    k = kappa(zu, zs)
    P = within(groups, rng, n)
    c = np.abs((zu[P] * np.conj(zs)[None]).sum(1)) / np.sqrt(np.sum(np.abs(zu) ** 2) * np.sum(np.abs(zs) ** 2))
    return float((1 + np.sum(c >= k)) / (1 + n))


def interval(est, draws):
    d = wrap(np.asarray(draws) - est) + est
    return float(np.quantile(d, .05)), float(np.quantile(d, .95))


def compartment(unit, wf, which, rng) -> dict:
    """All statistics of one compartment, waveform and unspliced definition."""
    obs, sp, un = unit["obs"], unit["spliced"], unit[which]
    grp, t = obs.group.to_numpy(), obs.t.to_numpy(float)
    L = projector(sp, wf)
    zs, zu = centred_states(sp, obs, wf, L), centred_states(un, obs, wf, L)
    row = dict(compartment=unit["name"], waveform=wf.name or "", unspliced=which, n=len(obs), primary=unit["primary"])
    # timed
    e = np.exp(-1j * W * t)
    bs, bu = np.sum(zs * e), np.sum(zu * e)
    T = t[within(grp, rng, N_PERM)]
    amp_null = np.abs((zs[None] * np.exp(-1j * W * T)).sum(1))
    row.update(rhythm_amp=abs(bs), rhythm_p=float((1 + np.sum(amp_null >= abs(bs))) / (1 + N_PERM)),
               delta_T=float(np.angle(bu * np.conj(bs)) * 24 / (2 * np.pi)))
    # untimed
    row.update(kappa=kappa(zu, zs), gate_p=gate_p(zu, zs, grp, rng, N_PERM), delta_U=lead(zu, zs))
    B = boot(grp, rng, N_BOOT)
    dU = np.array([lead(zu[b], zs[b]) for b in B])
    dT = np.array([np.angle(np.sum(zu[b] * e[b]) * np.conj(np.sum(zs[b] * e[b]))) * 24 / (2 * np.pi) for b in B])
    row["dU_lo"], row["dU_hi"] = interval(row["delta_U"], dU)
    row["dT_lo"], row["dT_hi"] = interval(row["delta_T"], dT)
    diff, diff_b = wrap(row["delta_U"] - row["delta_T"]), wrap(dU - dT)
    row.update(diff=float(diff), diff_lo=interval(diff, diff_b)[0], diff_hi=interval(diff, diff_b)[1])
    db = wrap(diff_b - diff) + diff
    row["tost_p"] = float(max(np.mean(db <= -1.5), np.mean(db >= 1.5)))
    logo = [lead(*(centred_states(x, obs, wf.drop([g]), projector(sp, wf.drop([g]))) for x in (un, sp)))
            for g in wf.genes]
    row["logo_min"], row["logo_max"] = float(np.min(logo)), float(np.max(logo))
    return row


def companions(unit, wf, which, seed) -> pd.DataFrame:
    obs, sp, un = unit["obs"], unit["spliced"], unit[which]
    grp = obs.group.to_numpy()
    lm = dev.log_mean(sp, sp.sum(1), min_count=0.0)  # pool: mean count >= 1, plus the waveform's own genes
    lm = lm[(sp.mean().reindex(lm.index) >= 1).to_numpy() | lm.index.isin(wf.genes)]
    rows = []
    for j, kind in enumerate(("matched", "responsive"), 1):
        for i, w in enumerate(pseudo_waveforms(lm, N_PSEUDO, np.random.default_rng([*seed, j]), wf, kind)):
            L = projector(sp, w)
            zs, zu = centred_states(sp, obs, w, L), centred_states(un, obs, w, L)
            rows.append(dict(kind=kind, set=i, kappa=kappa(zu, zs), delta_U=lead(zu, zs),
                             gate_p=gate_p(zu, zs, grp, np.random.default_rng([*seed, j, i]), N_PERM_PSEUDO)))
    return pd.DataFrame(rows)


def task(k):
    unit, wname, which = S["tasks"][k]
    wf = WAVEFORMS[wname]
    seed = [SEED, k]
    row = compartment(unit, wf, which, np.random.default_rng(seed))
    row["waveform"] = wname
    comp = companions(unit, wf, which, seed).assign(compartment=unit["name"], waveform=wname, unspliced=which)
    for kind in ("matched", "responsive"):
        c = comp[comp.kind == kind]
        ok = c[c.gate_p < .05]
        row[f"{kind}_kappa_q95"] = c.kappa.quantile(.95)
        row[f"{kind}_gate_share"] = len(ok) / len(c)
        row[f"{kind}_delta_q95"] = ok.delta_U.quantile(.95) if len(ok) else np.nan
    return row, comp


def gene_lags(unit) -> pd.DataFrame:
    """Descriptive: per clock gene, cosine-fit peak times of spliced and unspliced log fractions (group intercepts)."""
    obs, t = unit["obs"], unit["obs"].t.to_numpy(float)
    Z = np.column_stack([dev.design(obs, age=None, factors=("group",)), np.cos(W * t), np.sin(W * t)])
    M = np.column_stack([np.ones(len(t)), Z])
    out = []
    for g in WF.genes:
        if g not in unit["spliced"] or unit["unspliced"][g].mean() < 10:
            continue
        r = dict(compartment=unit["name"], gene=g, mean_unspliced=unit["unspliced"][g].mean())
        for k in ("spliced", "unspliced"):
            x = unit[k]
            y = np.log((x[g].to_numpy(float) + 0.5) / x.sum(1).to_numpy(float))
            b = np.linalg.lstsq(M, y, rcond=None)[0][-2:]
            r[f"peak_{k}"], r[f"amp_{k}"] = np.mod(np.arctan2(b[1], b[0]) / W, 24), np.hypot(*b)
        r["lag_h"] = float(wrap(r["peak_spliced"] - r["peak_unspliced"]))
        out.append(r)
    return pd.DataFrame(out)


def verdicts(t: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (wname, which), x in t[t.primary].groupby(["waveform", "unspliced"]):
        r = x[x.rhythm_p < .05].copy()
        if len(r):
            r["gate_holm"] = dev.holm(r.gate_p.to_numpy())
            r["H337_1"] = (r.gate_holm < .05) & (r.dU_lo > 0) & (r.dU_hi < 6)
            r["tost_holm"] = np.nan
            m = r.H337_1.to_numpy()
            if m.any():
                r.loc[m, "tost_holm"] = dev.holm(r.loc[m, "tost_p"].to_numpy())
            r["H337_2"] = r.H337_1 & (r.tost_holm < .05)
            r["H337_3"] = r.H337_1 & (r.delta_U > r.responsive_delta_q95)
        rows.append(pd.concat([r, x[x.rhythm_p >= .05].assign(H337_1=np.nan)]))
    v = pd.concat(rows)
    return v[["waveform", "unspliced", "compartment", "n", "rhythm_p", "kappa", "gate_p", "gate_holm", "delta_U",
              "dU_lo", "dU_hi", "delta_T", "dT_lo", "dT_hi", "diff", "diff_lo", "diff_hi", "tost_p", "tost_holm",
              "responsive_delta_q95", "H337_1", "H337_2", "H337_3"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sets", nargs="+", default=["baboon", "greenwell"])
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--reread", action="store_true", help="verdicts from the saved compartments.csv only")
    a = ap.parse_args()
    global N_PERM, N_PERM_PSEUDO, N_BOOT, N_PSEUDO
    limit_threads(1)
    t0, sfx = time.time(), "_".join(a.sets) + ("_check" if a.check else "")
    if a.reread:
        t = pd.read_csv(OUT / sfx / "compartments.csv")
        t["primary"] = ~t.compartment.isin(["baboon:Adrenal medulla", "greenwell:liver AR", "greenwell:liver DR"])
        t.to_csv(OUT / sfx / "compartments.csv", index=False)
        return report(t, OUT / sfx)
    units = [u for name in a.sets for u in load(name)]
    if a.check:  # no real statistic: unspliced := thinned spliced, times shuffled
        N_PERM, N_PERM_PSEUDO, N_BOOT, N_PSEUDO = 200, 50, 100, 10
        rng = np.random.default_rng(0)
        for u in units:
            u["unspliced"] = u["direct"] = u["spliced"].apply(lambda c: rng.binomial(c.astype(np.int64), 0.1))
            u["obs"] = u["obs"].assign(t=rng.permutation(u["obs"].t.to_numpy()))
    for u in units:
        print(u["name"], len(u["obs"]), "samples; groups", dict(u["obs"].group.value_counts()),
              "; clock genes present", sum(g in u["spliced"] for g in WF.genes), flush=True)
    S["tasks"] = [(u, w, which) for u in units for w in WAVEFORMS for which in ("unspliced", "direct")]
    out = OUT / sfx
    out.mkdir(parents=True, exist_ok=True)
    rows, comps = [], []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, (r, c) in enumerate(pool.imap_unordered(task, range(len(S["tasks"])))):
            rows.append(r)
            comps.append(c)
            print(f"{k + 1}/{len(S['tasks'])} {time.time() - t0:.0f}s", flush=True)
    t = pd.DataFrame(rows).sort_values(["waveform", "unspliced", "compartment"])
    t.to_csv(out / "compartments.csv", index=False)
    pd.concat(comps).to_csv(out / "companions.csv", index=False)
    pd.concat([gene_lags(u) for u in units]).to_csv(out / "genes.csv", index=False)
    report(t, out)
    print(f"runtime {time.time() - t0:.0f}s")


def report(t, out):
    """verdicts.csv and summary.csv: per waveform and intron definition, the shares of H337.1-3 over the rhythmic
    primary compartments and the expected lead (median and range of delta_T over them)."""
    v = verdicts(t)
    v.to_csv(out / "verdicts.csv", index=False)
    r = v[v.H337_1.notna()]
    s = r.groupby(["waveform", "unspliced"]).agg(
        rhythmic=("compartment", "size"), h1=("H337_1", "sum"), h2=("H337_2", "sum"), h3=("H337_3", "sum"),
        delta_T_median=("delta_T", "median"), delta_T_min=("delta_T", "min"), delta_T_max=("delta_T", "max")).reset_index()
    s.to_csv(out / "summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(v.round(3).to_string(index=False), s.round(3).to_string(index=False), sep="\n")


if __name__ == "__main__":
    main()
