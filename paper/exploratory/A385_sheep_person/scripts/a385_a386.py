"""A385 (animal loudness shared by organs, sheep) and A386 (glucocorticoid state against the animal's clock volume),
registered 2026-11-25.

SheepGTEx v0 (A263 genome-wide table), 35 peripheral tissues, adult and juvenile animals (metadata Table S1, corrected
stage); project = stratum (mixed strata). Per waveform (clock, 200 matched, 200 responsive; each on 14 genes and on the
core positions): one `Dataset` (every peripheral tissue an anchor); A385: 20 tissue halvings (`split_hosts`, seed
[385, 1]), leverage of the two half-host radii within each project with >= 30 animals, Fisher-z combined (weights n - 3),
median over halvings; A386: log radius of the full host per animal, regressed on base (number of tissues, project
indicators, cos / sin host phase), technical (median log depth, median mapping rate) and systemic (GC, IEG states) blocks.
Choice fixed in the power stage, before any value: the half-host radii are residualised on the number of the animal's
tissues in each half as well as on the phase harmonics (posterior shrinkage depends on it in both halves alike, which
would otherwise correlate the radii).

Stages:
  power  companions only (A385 leverages and A386 GC coefficients of the 400 sets, 14 genes and core positions), animal
         and project counts, analytic detection limit (A385) and power (A386). No clock value.
  check  2 sets, nothing written.
  score  A386 only (A385 not decidable): the clock fit on 14 genes and on the core positions; systemic and technical
         block partial R^2 against 10,000 permutations of the block rows (seed 386, as A360); standardised coefficients
         with animal-bootstrap 90 % intervals (2000); GC coefficient against the q05 of the stored companion
         coefficients (companions.csv, same model); demographic block (sex, stage) descriptive. Sensitivity
         (descriptive): the systemic permutation within project.

    RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A385_sheep_person/scripts/a385_a386.py --stage power --jobs 4
    RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A385_sheep_person/scripts/a385_a386.py --stage score
Outputs (`../output/`): companions.csv, power.csv; a386_tests.csv, a386_coefficients.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory" / "A263_companions" / "scripts"))
sys.path.insert(0, str(PAPER / "exploratory" / "A358_A360_gain_sources" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import score_a263 as A  # noqa: E402
from clocksync._core.controls import GC_TARGETS, IMMEDIATE_EARLY  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
META_XLSX = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/farm_gtex/sheep/raw/SheepGTEx_v0.MetaData.xlsx")
META = ["donor", "compartment", "depth", "stratum"]
WF = cs.gtex()
TRIO = ["PER1", "PER2", "NR1D1"]
CORE = WF.drop(TRIO)
CORE_IDX = [k for k, g in enumerate(WF.genes) if g in CORE.genes]
N_SETS, N_SPLITS, MIN_PROJECT, MIN_TISSUES = 200, 20, 30, 4
STAGES = ("Adult", "Juvenile")
S = {}


def metadata():
    m = pd.read_excel(META_XLSX, sheet_name="Table S1", header=1)
    m["Individual"] = m.Individual.astype(str)
    g = m.groupby("Individual")
    return pd.DataFrame({"stage": g["Developmental stage (corrected)"].first(), "sex": g["Sex (corrected)"].first(),
                         "mapping": pd.to_numeric(m["Mapping rate (%)"], errors="coerce").groupby(m.Individual).median()})


def load():
    pb = A.table("sheep")
    _, per, _ = A.groups("sheep", sorted(pb.compartment.unique()))
    md = metadata()
    keep = pb.compartment.isin(per) & pb.donor.astype(str).map(md.stage).isin(STAGES)
    pb = pb[keep].reset_index(drop=True)
    pb["donor"] = pb.donor.astype(str)
    return pb, per, md


def log_mean(pb):
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    return lm[~lm.index.duplicated()]


def leverage(a, b, cov):
    """corr of the two radii after residualising both on two phase harmonics of both states and on `cov` (n x k)."""
    ra, rb = np.abs(a), np.abs(b)
    H = np.column_stack([np.ones(len(a)), cov] + [g(k * p) for p in (np.angle(a), np.angle(b)) for k in (1, 2)
                                                  for g in (np.cos, np.sin)])
    res = lambda v: v - H @ np.linalg.lstsq(H, v, rcond=None)[0]  # noqa: E731
    return float(np.corrcoef(res(ra), res(rb))[0, 1])


def gc_states(pb, d, genes):
    """Per animal: mean over its tissues of the within-tissue standardised score (A359 rule, centred within tissue x
    project)."""
    key = pd.Series(np.arange(len(pb)), index=pd.MultiIndex.from_arrays([pb.donor, pb.compartment]))
    rows = key.loc[list(zip(d.resection.astype(str), d.compartment))].to_numpy()
    y = np.log((pb[genes].to_numpy(float)[rows] + 0.5) / pb.depth.to_numpy(float)[rows, None]).mean(1)
    grp = pd.Series(d.compartment.astype(str)) + "|" + pd.Series(d.group.astype(str))
    y = y - pd.Series(y).groupby(grp.to_numpy()).transform("mean").to_numpy()
    z = pd.Series(y).groupby(np.asarray(d.compartment)).transform(lambda v: (v - v.mean()) / v.std(ddof=1)).to_numpy()
    return pd.Series(z).groupby(np.asarray(d.resection).astype(str)).mean()


def design(h, animals):
    """A386 design (DataFrame) for the animals with >= MIN_TISSUES tissues: outcome log r and covariates."""
    a = animals.join(h[["x", "y", "radius"]], how="inner")
    a = a[a.n_tissues >= MIN_TISSUES].copy()
    a["log_r"] = np.log(a.radius)
    ph = np.angle(a.x + 1j * a.y)
    a["ph_cos"], a["ph_sin"] = np.cos(ph), np.sin(ph)
    return a


def gc_coef(a, projects):
    cols = ["n_tissues", "ph_cos", "ph_sin", "log_depth", "mapping", "gc", "ieg"]
    a = a.dropna(subset=["log_r"] + cols)
    Z = (a[cols] - a[cols].mean()) / a[cols].std(ddof=1)
    P = pd.get_dummies(a.project, drop_first=True).to_numpy(float)
    y = (a.log_r - a.log_r.mean()) / a.log_r.std(ddof=1)
    X = np.column_stack([np.ones(len(a)), P, Z.to_numpy()])
    beta = np.linalg.lstsq(X, y.to_numpy(), rcond=None)[0]
    k = X.shape[1] - 2  # position of the gc term
    se = float(np.sqrt(np.linalg.inv(X.T @ X)[k, k]))  # outcome variance 1 (power planning)
    return float(beta[k]), se, len(a)


BLOCKS = {"technical": ["log_depth", "mapping"], "systemic": ["gc", "ieg"], "demographic": ["male", "adult"]}


def _X(a, cols):
    """Intercept, project indicators and the standardised `cols`."""
    Z = ((a[cols] - a[cols].mean()) / a[cols].std(ddof=1)).to_numpy(float)
    P = pd.get_dummies(a.project, drop_first=True).to_numpy(float)
    return np.column_stack([np.ones(len(a)), P, Z]), P.shape[1] + 1


def _rss(y, X):
    return float(((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2).sum())


def partial_r2(a, block, given, rng, n=10_000, within=False):
    """Partial R^2 of `block` given project + `given`; p from n permutations of the block rows (within project if
    `within`)."""
    y = a.log_r.to_numpy(float)
    Xg, _ = _X(a, given)
    Xb = ((a[block] - a[block].mean()) / a[block].std(ddof=1)).to_numpy(float)
    r_red = _rss(y, Xg)
    pr2 = 1 - _rss(y, np.column_stack([Xg, Xb])) / r_red
    groups = [np.flatnonzero(a.project.to_numpy() == p) for p in a.project.unique()]
    null = np.empty(n)
    for k in range(n):
        if within:
            idx = np.arange(len(y))
            for g in groups:
                idx[g] = g[rng.permutation(len(g))]
        else:
            idx = rng.permutation(len(y))
        null[k] = 1 - _rss(y, np.column_stack([Xg, Xb[idx]])) / r_red
    return pr2, float((1 + np.sum(null >= pr2)) / (n + 1))


def std_coefs(a, cols):
    X, k0 = _X(a, cols)
    y = ((a.log_r - a.log_r.mean()) / a.log_r.std(ddof=1)).to_numpy(float)
    return np.linalg.lstsq(X, y, rcond=None)[0][k0:]


def score(pbk, per, animals, gc, ieg, project):
    """A386 scoring: clock on 14 genes and on the core positions."""
    C = pd.read_csv(OUT / "companions.csv")
    base = ["n_tissues", "ph_cos", "ph_sin"]
    full = base + BLOCKS["technical"] + BLOCKS["systemic"]
    tests, coefs = [], []
    for gs, w, hyp in (("14 genes", WF, "H386.1"), ("core", CORE, "H386.2")):
        rng = np.random.default_rng(386)
        ds = cs.Dataset(pbk[META + list(w.genes)], anchors=per, waveform=w, mixed_strata=True, name=f"a386 {gs}",
                        min_donors=5)
        d = ds.data
        tis = pd.Series(np.asarray(d.resection).astype(str))
        h = ds.host(per)
        h.index = h.index.astype(str)
        an = animals.assign(n_tissues=tis.value_counts())
        an["gc"], an["ieg"] = gc_states(pbk, d, gc), gc_states(pbk, d, ieg)
        a = design(h, an).dropna(subset=["log_r"] + full)
        beta, _, n = gc_coef(a, project)
        r2_sys, p_sys = partial_r2(a, BLOCKS["systemic"], base + BLOCKS["technical"], rng)
        _, p_sys_w = partial_r2(a, BLOCKS["systemic"], base + BLOCKS["technical"], rng, within=True)
        r2_tech, p_tech = partial_r2(a, BLOCKS["technical"], base + BLOCKS["systemic"], rng)
        c = C[C.label == gs]
        q05 = {k: float(np.nanquantile(c[c.kind == k].gc_coef, 0.05)) for k in ("matched", "responsive")}
        rank = {k: float((c[c.kind == k].gc_coef <= beta).mean()) for k in ("matched", "responsive")}
        tests.append(dict(test=hyp, genes=gs, animals=n, projects=a.project.nunique(), gc_coef=beta,
                          systemic_partial_r2=r2_sys, p_systemic=p_sys, p_systemic_within_project=p_sys_w,
                          technical_partial_r2=r2_tech, p_technical=p_tech, matched_q05=q05["matched"],
                          responsive_q05=q05["responsive"], frac_matched_le=rank["matched"],
                          frac_responsive_le=rank["responsive"],
                          passes=bool(p_sys < 0.05 and beta < 0 and beta < q05["matched"] and beta < q05["responsive"])))
        b = std_coefs(a, full)
        boot = []
        for _ in range(2000):
            s = a.iloc[rng.integers(0, len(a), len(a))]
            try:
                boot.append(std_coefs(s, full))
            except (ValueError, np.linalg.LinAlgError):
                continue
        boot = np.array(boot)
        coefs.append(pd.DataFrame({"genes": gs, "model": "registered", "term": full, "std_coef": b,
                                   "lo90": np.nanquantile(boot, 0.05, 0), "hi90": np.nanquantile(boot, 0.95, 0),
                                   "n_boot": len(boot)}))
        ad = a.assign(male=a.sex.map({"Male": 1.0, "Female": 0.0}),
                      adult=(a.stage == "Adult").astype(float)).dropna(subset=BLOCKS["demographic"])
        r2_dem, p_dem = partial_r2(ad, BLOCKS["demographic"], full, rng)
        tests.append(dict(test="demographic (descriptive)", genes=gs, animals=len(ad), systemic_partial_r2=r2_dem,
                          p_systemic=p_dem))
        coefs.append(pd.DataFrame({"genes": gs, "model": "with demographic (descriptive)",
                                   "term": full + BLOCKS["demographic"],
                                   "std_coef": std_coefs(ad, full + BLOCKS["demographic"])}))
        print(f"{gs}: {n} animals done", flush=True)
    T, K = pd.DataFrame(tests), pd.concat(coefs, ignore_index=True)
    T.to_csv(OUT / "a386_tests.csv", index=False)
    K.to_csv(OUT / "a386_coefficients.csv", index=False)
    pd.set_option("display.width", 250)
    print(T.round(4).to_string(index=False))
    print(K.round(3).to_string(index=False))


def halvings(per):
    """The 20 tissue halvings (as `elasticity.split_hosts`, seed [385, 1])."""
    rng = np.random.default_rng([385, 1])
    out = []
    for _ in range(N_SPLITS):
        p = list(rng.permutation(per))
        out.append((sorted(p[: len(p) // 2]), sorted(p[len(p) // 2:])))
    return out


def project_leverage(ds, tis, halves, project, rng=None):
    """Median over halvings of the Fisher-z combined within-project leverage (projects with >= MIN_PROJECT animals in
    both halves) and the median total weight sum(n - 3). rng: permute the animals of half B within project (null)."""
    L, N = [], []
    for sa, sb in halves:
        hA, hB = ds.host(sa), ds.host(sb)
        j = hA.index.astype(str).intersection(hB.index.astype(str))
        hA.index, hB.index = hA.index.astype(str), hB.index.astype(str)
        a, b = (hA.x + 1j * hA.y)[j].to_numpy(), (hB.x + 1j * hB.y)[j].to_numpy()
        cnt = [tis[tis.comp.isin(s)].groupby("donor").size().reindex(j).to_numpy(float) for s in (sa, sb)]
        cov = np.log(np.column_stack(cnt))
        pr = project.reindex(j).to_numpy()
        z = w = 0.0
        for p in np.unique(pr):
            m = np.flatnonzero(pr == p)
            if len(m) < MIN_PROJECT:
                continue
            mb = m if rng is None else rng.permutation(m)
            z += (len(m) - 3) * np.arctanh(leverage(a[m], b[mb], cov[m]))
            w += len(m) - 3
        L.append(np.tanh(z / w) if w > 0 else np.nan)
        N.append(w)
    return float(np.nanmedian(L)), float(np.median(N))


def task(i):
    warnings.filterwarnings("ignore")
    lab, w = S["wfs"][i]
    pb, per = S["pb"], S["per"]
    ds = cs.Dataset(pb[META + list(w.genes)], anchors=per, waveform=w, mixed_strata=True, name=str(lab), min_donors=5)
    d = ds.data
    tis = pd.DataFrame({"donor": np.asarray(d.resection).astype(str), "comp": np.asarray(d.compartment)})
    lev, wn = project_leverage(ds, tis, halvings(per), S["project"])
    h = ds.host(per)
    h.index = h.index.astype(str)
    animals = S["animals"].assign(n_tissues=tis.donor.value_counts())
    animals["gc"], animals["ieg"] = gc_states(pb, d, S["gc"]), gc_states(pb, d, S["ieg"])
    beta, se, n = gc_coef(design(h, animals), S["project"])
    return dict(id=i, label=str(lab[0]), kind=lab[1], set=lab[2], leverage=lev, weight_n3=wn, gc_coef=beta,
                gc_se_design=se, animals_a386=n)


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "check", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    pb, per, md = load()
    lm = log_mean(pb)
    n = 2 if a.stage == "check" else N_SETS
    sets = {kind: pseudo_waveforms(lm, n, np.random.default_rng([385, tag]), WF, kind=kind)
            for kind, tag in (("matched", 2), ("responsive", 3))}
    core_of = lambda w, j: Waveform(tuple(w.genes[k] for k in CORE_IDX), CORE.amp, CORE.psi, None, None, f"core{j}")  # noqa: E731
    wfs = [((gs, kind, j), w if gs == "14 genes" else core_of(w, j))
           for kind, ws in sets.items() for j, w in enumerate(ws) for gs in ("14 genes", "core")]
    clock = {canonical(g) for g in WF.genes}
    gc = [g for g in dict.fromkeys(GC_TARGETS) if g in pb.columns and g not in {"NFIL3", "BHLHE40"} and canonical(g) not in clock]
    ieg = [g for g in dict.fromkeys(IMMEDIATE_EARLY) if g in pb.columns and canonical(g) not in clock]
    keep = sorted({g for _, w in wfs for g in w.genes} | set(WF.genes) | set(gc) | set(ieg))
    pbk = pb[META + keep].copy()
    project = pbk.groupby("donor").stratum.first()
    dep = np.log(pbk.depth).groupby(pbk.donor).median()
    animals = pd.DataFrame({"project": project, "log_depth": dep}).join(md[["mapping", "stage", "sex"]])
    S.update(wfs=wfs, pb=pbk, per=per, project=project, animals=animals, gc=gc, ieg=ieg)
    del pb
    if a.stage == "score":
        score(pbk, per, animals, gc, ieg, project)
        print(f"({time.time() - t0:.0f}s)")
        return
    print(f"{len(per)} tissues, {pbk.donor.nunique()} animals, {project.nunique()} projects, GC {len(gc)} IEG {len(ieg)} "
          f"genes, {len(wfs)} waveforms ({time.time() - t0:.0f}s)", flush=True)
    rows = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, r in enumerate(pool.imap_unordered(task, range(len(wfs)))):
            rows.append(r)
            if (k + 1) % 25 == 0 or k + 1 == len(wfs):
                print(f"{k + 1}/{len(wfs)} ({time.time() - t0:.0f}s)", flush=True)
    C = pd.DataFrame(rows).sort_values("id")
    if a.stage == "check":
        print(C[["label", "kind", "weight_n3", "animals_a386"]].to_string(), f"({time.time() - t0:.0f}s)")
        return
    OUT.mkdir(parents=True, exist_ok=True)
    C.to_csv(OUT / "companions.csv", index=False)
    rows = []
    for gs in ("14 genes", "core"):
        c = C[C.label == gs]
        q95 = {k: float(np.nanquantile(c[c.kind == k].leverage, 0.95)) for k in ("matched", "responsive")}
        se_z = 1 / np.sqrt(float(c.weight_n3.median()))
        thr = max(q95.values())
        l80 = float(np.tanh(np.arctanh(thr) + (1.645 + 0.84) * se_z))
        q05 = {k: float(np.nanquantile(c[c.kind == k].gc_coef, 0.05)) for k in ("matched", "responsive")}
        se_gc = float(c.gc_se_design.median())
        g_thr = min(q05["matched"], q05["responsive"], -1.96 * se_gc)
        pw = lambda b: float(1 - __import__("scipy").stats.norm.cdf((b - g_thr) / se_gc))  # noqa: E731
        rows.append(dict(genes=gs, weight_n3=float(c.weight_n3.median()), animals_a386=int(c.animals_a386.median()),
                         lev_matched_median=float(c[c.kind == "matched"].leverage.median()),
                         lev_matched_q95=q95["matched"], lev_responsive_q95=q95["responsive"], L80=l80,
                         gate_A385=bool(l80 <= 0.40),
                         gc_matched_median=float(c[c.kind == "matched"].gc_coef.median()), gc_matched_q05=q05["matched"],
                         gc_responsive_q05=q05["responsive"], gc_se=se_gc, gc_threshold=g_thr,
                         power_gc_0235=pw(-0.235), power_gc_0118=pw(-0.1175), gate_A386=bool(pw(-0.235) >= 0.8),
                         failed=int(C[C.label == gs].leverage.isna().sum())))
    P = pd.DataFrame(rows)
    P.to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 250)
    print(P.round(3).to_string(index=False), f"\n({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
