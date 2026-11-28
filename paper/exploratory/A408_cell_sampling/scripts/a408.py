"""A408 (is a class's own clock departure the sampling of its incoherent cells?), registered 2026-11-27.

HLCA healthy lung (A17 `lung_hlca`, learned anchors), cells of the A17 healthy donors under the A17 lung class rule and
majority-assay filter (`lung.py`, `census.py`; the per-cell clock-gene file `$SCCC_DATA/cells/hlca_full`). Per target
class x donor: 20 draws of n in {25, 50, 100, 200, 400} cells without replacement, pseudobulked and read in place of the
target compartment (host = full pseudobulks of the other anchors, `Dataset.fit`); read-thinning control (the 400-cell draw
binomially thinned to the depth of the n-cell draw); Delta(n) = psi_cells(n) - psi_reads(n); OLS of Delta on 1/n through the
origin (s2); f = (s2 / 400) / psi_cells(400); 50 expression-matched companion sets (A263 rule, seed [408, 1]); random-effects
pooling over units (`clocksync.stats.meta`), donor bootstraps. Notes 1 and 2 of the registration list the choices.

Stages:
  prep   units (class x donor counts under the A17 filters), the 50 companion sets, the cell request for the Census query
         (`a408_census.py`, run in `census_env`). No clock psi.
  power  unit table, companion psi at 400 cells (never the real clock psi), and the simulation calibration: clock counts of
         the real cells replaced by Poisson counts from phase-scattered cells (cell part of the own term at 400 cells = 0, 0.5
         and 1 x the companions' psi_cells(400)), read through the same pipeline: bias of f and donor-bootstrap SE.
  score  clock and companions through the pipeline; H408.1, H408.2 (implemented, run only after the power note).

    paper/common/run.sh paper/exploratory/A408_cell_sampling/scripts/a408.py --stage power --jobs 4
Outputs (`../output/`): units.csv, sets.csv, power.csv, sim.csv, companions_c400.csv; tests.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "common" / "datasets" / "_atlas"))
sys.path.insert(0, str(PAPER / "results" / "figure_2" / "scripts" / "atlas_fits"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.anchors import family  # noqa: E402
from common.anchors import anchors as a17_anchors  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 408
NS = (25, 50, 100, 200, 400)
N_DRAWS, MIN_CELLS, MIN_DONORS = 20, 400, 8
TARGET_FAMILIES = {"fibroblast_stromal", "endothelial", "mural", "epithelial"}  # parenchyma: lung epithelium
WF = cs.gtex()
CLOCK = list(WF.genes)
CELLS = DATA / "cells" / "hlca_full" / "cells.parquet"
COMP_DIR = DATA / "cells" / "hlca_a408_companions"
REQ = OUT / "census_request"
ATLAS = "lung_hlca"


# ---- cells under the A17 rules -------------------------------------------------------------------------------------

def cells_a17(columns=("observation_joinid", "donor_id", "assay", "ann_level_2", "ann_level_3", "tissue",
                       "tissue_sampling_method", "cell_type", "study", "raw_sum", *CLOCK)) -> pd.DataFrame:
    """Per-cell table of the cells that enter the A17 `lung_hlca` table: lung sampling methods (no nose), class from
    `ann_level_3` (smooth muscle from `ann_level_2`), healthy donors of the fit, the donor's majority assay."""
    import lung as LU
    c = pd.read_parquet(CELLS, columns=list(columns))
    # the group pseudobulk of the atlas (`atlas/hlca_full`) keeps groups of >= 10 cells (donor x tissue x cell type x assay x
    # study x sampling method x ann_level_2 x ann_level_3; `tissue_general` is constant within tissue)
    gk = ["donor_id", "tissue", "cell_type", "assay", "study", "tissue_sampling_method", "ann_level_2", "ann_level_3"]
    c = c[c.groupby(gk, observed=True).raw_sum.transform("size") >= 10]
    cls = c.ann_level_3.astype(str).map(LU.L3)
    cls[c.ann_level_2.astype(str) == "Smooth muscle"] = "smooth_muscle"
    keep = (c.tissue_sampling_method.astype(str).isin(LU.HLCA_TISSUE) & (c.tissue.astype(str) != "nose")
            & cls.notna()).to_numpy()
    c = c[keep].assign(compartment="lung:" + cls[keep].astype(str), donor="hlca_full:" + c.donor_id[keep].astype(str))
    # donor set of the A17 healthy fit
    dn = pd.read_parquet(DATA / "pseudobulk" / "lung_hlca_donors.parquet")
    ld = LU.hlca_donors()
    ok = ld.index[~ld.lung_condition.str.contains(";")]
    c = c[c.donor.isin(dn.index[dn.healthy.astype(bool)])]
    top = c.groupby(["donor", "assay"]).size().reset_index(name="n").sort_values(["donor", "n", "assay"])
    maj = top.drop_duplicates("donor", keep="last").set_index("donor").assay
    c = c[c.assay.to_numpy() == maj.reindex(c.donor).to_numpy()]
    return c.drop(columns=["donor_id", "ann_level_2", "ann_level_3", "tissue", "tissue_sampling_method", "cell_type", "study"]).reset_index(drop=True)


def a17_table():
    """The A17 `lung_hlca` healthy table (`Target`) and its anchors (learned rule)."""
    from _shared import Target
    T = Target(ATLAS)
    anc = [a for a in a17_anchors(ATLAS, "learned") if a in set(T.pb.compartment)]
    return T, anc


def units(c: pd.DataFrame, T, anc) -> pd.DataFrame:
    """Class x atlas units: target family, >= 8 donors with >= 400 cells, host >= 3 other anchors."""
    n = c.groupby(["compartment", "donor"]).size().rename("cells").reset_index()
    fam = {k: family(k) for k in n.compartment.unique()}
    rows = []
    for k, g in n.groupby("compartment"):
        d400 = int((g.cells >= MIN_CELLS).sum())
        host = len([a for a in anc if a != k])
        rows.append(dict(atlas=ATLAS, compartment=k, family=fam[k], donors_cells=len(g), donors_400=d400,
                         median_cells=float(g.cells.median()), host_others=host,
                         target=fam[k] in TARGET_FAMILIES, anchor=k in anc,
                         enters=bool(fam[k] in TARGET_FAMILIES and d400 >= MIN_DONORS and host >= 3)))
    return pd.DataFrame(rows)


def check_against_pb(c: pd.DataFrame, T, n=300, seed=0) -> dict:
    """Cell-level sums against the A17 table rows (counts of the clock genes, n_cells, depth)."""
    rng = np.random.default_rng(seed)
    s = c.groupby(["donor", "compartment"]).agg(n_cells=("raw_sum", "size"), depth=("raw_sum", "sum"),
                                                **{g: (g, "sum") for g in CLOCK})
    pb = T.pb.set_index(["donor", "compartment"])
    both = s.index.intersection(pb.index)
    both = both[rng.permutation(len(both))[:n]]
    a, b = s.loc[both], pb.loc[both]
    return dict(rows_cells=len(s), rows_pb=len(pb), common=len(s.index.intersection(pb.index)), checked=len(both),
                max_abs_clock=float(np.abs(a[CLOCK].to_numpy(float) - b[CLOCK].to_numpy(float)).max()),
                max_abs_ncells=float(np.abs(a.n_cells - b.n_cells).max()),
                max_rel_depth=float((np.abs(a.depth - b.depth) / b.depth).max()))


# ---- companions -----------------------------------------------------------------------------------------------------

def genome_table():
    """Genome-wide donor x compartment table of the healthy HLCA (`lung_hlca_genome`, the rows of the A17 table)."""
    from common import data
    pb, _ = data._builder("lung_hlca_genome").build()
    return pb


def log_mean(pb: pd.DataFrame) -> pd.Series:
    """A263 matching variable: log mean fraction over all rows of the genome table, genes with mean count >= 1."""
    genes = [c for c in pb.columns if c not in ("donor", "compartment", "stratum", "n_cells", "depth")]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)
    return lm[x.mean(0) >= 1]


def companion_sets(lm: pd.Series, n_sets: int = 50):
    from clocksync.paired import pseudo_waveforms
    return pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, 1]), WF, "matched")


# ---- cell pools, draws, pseudobulks ---------------------------------------------------------------------------------

class Pool:
    """Cells of one unit (class x atlas): donors with >= 400 cells, per-cell depth, clock counts and (when loaded) the
    companion genes' counts (sparse, columns `comp_genes`)."""

    def __init__(self, comp, ci, cells: pd.DataFrame, T, comp_counts=None):
        d = cells[cells.compartment == comp]
        n = d.groupby("donor").size()
        d = d[d.donor.isin(n[n >= MIN_CELLS].index)].sort_values(["donor", "observation_joinid"]).reset_index(drop=True)
        self.comp, self.ci = comp, ci
        self.donors = sorted(d.donor.unique())
        self.code = pd.Categorical(d.donor, categories=self.donors).codes
        self.start = np.searchsorted(self.code, np.arange(len(self.donors) + 1))  # cells of donor j: start[j]:start[j+1]
        self.raw = d.raw_sum.to_numpy(float)
        self.clock = d[CLOCK].to_numpy(np.int32)
        self.ids = d.observation_joinid.to_numpy()
        row = T.pb[T.pb.compartment == comp].set_index("donor")
        self.stratum = row.stratum.reindex(self.donors).astype(str).to_numpy()
        self.comp_x = None if comp_counts is None else comp_counts(self.ids)

    def draw(self, n: int, r: int):
        """Cell indices (concatenated over donors, donor code per index) of draw r of size n."""
        rng = np.random.default_rng([SEED, 0, self.ci, n, r])
        idx = [self.start[j] + rng.choice(self.start[j + 1] - self.start[j], n, replace=False)
               for j in range(len(self.donors))]
        return np.concatenate(idx), np.repeat(np.arange(len(self.donors)), n)

    def pseudobulk(self, X, n: int, r: int):
        """Summed counts (donors x genes) and depth of draw (n, r); X: cells x genes (dense or CSR)."""
        i, g = self.draw(n, r)
        M = sp.csr_matrix((np.ones(len(i)), (g, np.arange(len(i)))), shape=(len(self.donors), len(i)))
        out = M @ X[i]
        return (out.toarray() if sp.issparse(out) else np.asarray(out)).astype(float), M @ self.raw[i]


def thin_reads(counts, depth400, depth_n, seed):
    """Binomial thinning of the 400-cell pseudobulk to the depth of the n-cell draw (p = depth_n / depth_400, at most 1)."""
    rng = np.random.default_rng(seed)
    p = np.minimum(depth_n / depth400, 1.0)
    return rng.binomial(counts.astype(np.int64), p[:, None]).astype(float), p * depth400


def target_tables(pool: Pool, X, genes):
    """Cell-draw and read-thinned pseudobulks of every (n, r): {(kind, n, r): (counts DataFrame, depth)}; X holds `genes`."""
    out, base = {}, {}
    for r in range(N_DRAWS):
        base[r] = pool.pseudobulk(X, 400, r)
    for n in NS:
        for r in range(N_DRAWS):
            c, d = pool.pseudobulk(X, n, r) if n < 400 else base[r]
            out[("cells", n, r)] = (c, d)
            if n < 400:
                out[("reads", n, r)] = thin_reads(base[r][0], base[r][1], d, [SEED, 0, pool.ci, n, r, 1])
            else:
                out[("reads", n, r)] = base[r]
    return out


def build_dataset(host_rows: pd.DataFrame, pool: Pool, tabs, genes, waveform, T, anc_k):
    """Host rows (donor, compartment, stratum, depth, genes) + one target compartment per (kind, n, r)."""
    parts = [host_rows[["donor", "compartment", "stratum", "depth", *genes]]]
    for (kind, n, r), (c, d) in tabs.items():
        t = pd.DataFrame(c, columns=genes)
        t.insert(0, "depth", d)
        t.insert(0, "stratum", pool.stratum)
        t.insert(0, "compartment", f"{kind}|{n}|{r}")
        t.insert(0, "donor", pool.donors)
        parts.append(t)
    return cs.Dataset(pd.concat(parts, ignore_index=True), anchors=anc_k, donors=T.donors, waveform=waveform,
                      name=f"a408|{pool.comp}")


def psi_of(ds, names) -> pd.Series:
    """psi of the target compartments `names` on the host of the anchors (`Dataset.fit`, 14 genes, no bootstrap)."""
    try:
        return ds.fit(n_boot=0, gate=False, targets=list(names)).table.psi.reindex(names)
    except (ValueError, np.linalg.LinAlgError):
        return pd.Series(np.nan, index=names)


def estimands(psi: pd.Series) -> dict:
    """psi (index 'kind|n|r') -> Delta(n), s2 (OLS of Delta on 1/n through the origin), f; s2_alt regresses on 1/n - 1/400."""
    P = psi.rename_axis("k").reset_index()
    P[["kind", "n", "r"]] = P.k.str.split("|", expand=True)
    P["n"], P["r"] = P.n.astype(int), P.r.astype(int)
    w = P.pivot_table(index=["n", "r"], columns="kind", values="psi")
    delta = (w["cells"] - w["reads"]).groupby("n").mean()
    n = delta.index.to_numpy(float)
    x = 1.0 / n
    s2 = float(np.sum(x * delta.to_numpy()) / np.sum(x * x))
    m = n < 400
    xa = x[m] - 1 / 400
    s2a = float(np.sum(xa * delta.to_numpy()[m]) / np.sum(xa * xa))
    p400 = float(w.loc[400, "cells"].mean())
    p200 = float(w.loc[200, "cells"].mean())
    return dict(s2=s2, f=s2 / 400 / p400 if p400 > 0 else np.nan, s2_alt=s2a,
                f_alt=s2a / 400 / p400 if p400 > 0 else np.nan, psi400=p400, psi200=p200,
                **{f"delta{int(k)}": float(v) for k, v in delta.items()})


# ---- shared state of the workers ------------------------------------------------------------------------------------
G = {}


def init_worker(state):
    G.update(state)
    limit_threads(1)
    warnings.filterwarnings("ignore")


def comp_loader():
    """ids -> CSR of the companion genes' counts (cells of the Census query), and the gene order."""
    cells = pd.read_csv(COMP_DIR / "cells.csv").observation_joinid.astype(str)
    X = sp.load_npz(COMP_DIR / "counts.npz").tocsr()
    genes = pd.read_csv(COMP_DIR / "genes.csv", header=None)[0].astype(str).tolist()
    pos = pd.Series(np.arange(len(cells)), index=cells.to_numpy())
    return (lambda ids: X[pos.reindex(ids).to_numpy()]), genes


def host_rows_of(table: pd.DataFrame, anc_k, genes):
    return table[table.compartment.isin(anc_k)][["donor", "compartment", "stratum", "depth", *genes]]


def setup(need_comp=True):
    """Cells, A17 table, units, pools (clock and, when needed, companion counts), companion sets and genome host rows."""
    t0 = time.time()
    c = cells_a17()
    T, anc = a17_table()
    U = units(c, T, anc)
    ent = U[U.enters].compartment.tolist()
    st = dict(T=T, anc=anc, U=U, ent=ent)
    if need_comp:
        load, cgenes = comp_loader()
        st["pools"] = {k: Pool(k, i, c, T, load) for i, k in enumerate(ent)}
        st["cgenes"] = cgenes
        gt = genome_table()
        st["sets"] = companion_sets(log_mean(gt))
        gt = gt[gt.compartment.isin(anc)][["donor", "compartment", "stratum", "depth", *cgenes]].copy()
        # one stratum per donor, that of the A17 table (the genome table writes the same assay | study label differently)
        sd = T.pb.groupby("donor").stratum.first()
        gt = gt[gt.donor.isin(sd.index)]
        gt["stratum"] = gt.donor.map(sd).to_numpy()
        st["gt"] = gt
    else:
        st["pools"] = {k: Pool(k, i, c, T) for i, k in enumerate(ent)}
    print(f"setup {time.time() - t0:.0f}s: units {ent}", flush=True)
    return st


def anchors_without(k):
    return [a for a in G["anc"] if a != k]


# ---- tasks ------------------------------------------------------------------------------------------------------------

def task_companion(job):
    """psi of the 20 400-cell draws of one unit for one companion set: mean and n valid (no n < 400, no reads)."""
    k, j = job
    pool, w = G["pools"][k], G["sets"][j]
    cols = [G["cgenes"].index(g) for g in w.genes]
    tabs = {("cells", 400, r): (c[:, cols], d) for r in range(N_DRAWS)
            for c, d in [pool.pseudobulk(G["pools"][k].comp_x, 400, r)]}
    ds = build_dataset(host_rows_of(G["gt"], anchors_without(k), w.genes), pool, tabs, list(w.genes), w, G["T"],
                       anchors_without(k))
    psi = psi_of(ds, [f"cells|400|{r}" for r in range(N_DRAWS)])
    return dict(compartment=k, set=j, psi400=float(psi.mean()), valid=int(psi.notna().sum()))


def host_states(k):
    """Host state (x, y) per donor of the anchors without k (host of the other compartments only)."""
    T, anc = G["T"], anchors_without(k)
    tab = T.pb[T.pb.compartment.isin(anc + ["lung:T"])]
    ds = cs.Dataset(tab, anchors=anc, donors=T.donors, name="host")
    return ds.fit(n_boot=0, gate=False, targets=["lung:T"]).host[["x", "y"]]


def simulate_cells(pool: Pool, h: pd.DataFrame, P: float, kfac: float, rng):
    """Per-cell clock counts: Poisson(depth_c * p_g * exp(X_g . beta_c)), beta_c = host state of the donor + donor-level own
    term (variance psi_inf per dimension) + the cell's incoherent part (uniform phase, amplitude sqrt(2 s2): variance s2 per
    dimension). s2 = kfac * P (P: companions' psi_cells(400)); psi_inf = P - s2 / 400."""
    s2 = kfac * P
    psi_inf = max(P - s2 / 400, 0.0)
    X = WF.X
    frac = pool.clock.sum(0) / pool.raw.sum()
    H = h.reindex(pool.donors).fillna(0.0).to_numpy(float)  # donors without a host state are not read by the fit
    u = rng.normal(0, np.sqrt(psi_inf), (len(pool.donors), 2))
    th = rng.uniform(0, 2 * np.pi, len(pool.raw))
    e = np.sqrt(2 * s2) * np.column_stack([np.cos(th), np.sin(th)])
    beta = (H + u)[pool.code] + e
    return rng.poisson(pool.raw[:, None] * frac[None, :] * np.exp(beta @ X.T)).astype(np.int32), s2, psi_inf


def run_pipeline(pool: Pool, X, genes, host_rows, waveform, anc_k):
    tabs = target_tables(pool, X, genes)
    ds = build_dataset(host_rows, pool, tabs, genes, waveform, G["T"], anc_k)
    names = [f"{kind}|{n}|{r}" for (kind, n, r) in tabs]
    return ds, names


def boot_f(ds, names, B, rng):
    out = []
    for _ in range(B):
        out.append(estimands(psi_of(ds.resample(rng), names))["f"])
    return np.array(out)


def task_sim(job):
    """One simulation replicate: unit k, cell-part level kfac (0, 0.5, 1), replicate rep; f and its donor-bootstrap SE."""
    k, kfac, rep, P, B = job
    pool = G["pools"][k]
    anc_k = anchors_without(k)
    rng = np.random.default_rng([SEED, 9, pool.ci, int(round(kfac * 10)), rep])
    h = host_states(k)
    Xs, s2, psi_inf = simulate_cells(pool, h, P, kfac, rng)
    t0 = time.time()
    ds, names = run_pipeline(pool, Xs, CLOCK, host_rows_of(G["T"].pb, anc_k, CLOCK), WF, anc_k)
    est = estimands(psi_of(ds, names))
    t1 = time.time()
    fb = boot_f(ds, names, B, np.random.default_rng([SEED, 10, pool.ci, int(round(kfac * 10)), rep])) if B else np.array([])
    f_true = s2 / 400 / (psi_inf + s2 / 400)
    return dict(compartment=k, kfac=kfac, rep=rep, P=P, s2_true=s2, psi_inf_true=psi_inf, f_true=f_true, **est,
                f_se=float(np.nanstd(fb, ddof=1)) if len(fb) > 2 else np.nan, boot_valid=int(np.isfinite(fb).sum()),
                secs_point=t1 - t0, secs_boot=time.time() - t1)


# ---- stages ---------------------------------------------------------------------------------------------------------

def stage_prep(a):
    OUT.mkdir(exist_ok=True)
    t0 = time.time()
    c = cells_a17()
    T, anc = a17_table()
    chk = check_against_pb(c, T)
    print("check vs A17 table:", chk, f"{time.time() - t0:.0f}s", flush=True)
    U = units(c, T, anc)
    U.to_csv(OUT / "units.csv", index=False)
    print(U.round(1).to_string(index=False), flush=True)
    pd.DataFrame([chk]).to_csv(OUT / "check_cells_vs_a17.csv", index=False)
    gt = genome_table()
    lm = log_mean(gt)
    sets = companion_sets(lm)
    pd.DataFrame([dict(set=i, **{f"g{j}": g for j, g in enumerate(w.genes)}) for i, w in enumerate(sets)]).to_csv(
        OUT / "sets.csv", index=False)
    genes = sorted({g for w in sets for g in w.genes})
    ent = U[U.enters]
    cc = c[c.compartment.isin(ent.compartment)]
    d400 = cc.groupby(["compartment", "donor"]).size()
    d400 = d400[d400 >= MIN_CELLS].reset_index()[["compartment", "donor"]]
    cc = cc.merge(d400, on=["compartment", "donor"])
    REQ.mkdir(exist_ok=True)
    cc[["observation_joinid"]].to_csv(REQ / "cells.csv", index=False)
    pd.Series(genes).to_csv(REQ / "genes.csv", index=False, header=False)
    print(f"sets 50, genes {len(genes)}, cells requested {len(cc)} of {ent.compartment.nunique()} units "
          f"{time.time() - t0:.0f}s", flush=True)


def pool_units(f: pd.DataFrame, se: pd.DataFrame | pd.Series):
    """Random-effects pooling (`clocksync.stats.meta`) of one value per unit with its SE."""
    ok = np.isfinite(f) & np.isfinite(se) & (se > 0)
    if ok.sum() < 2:
        return dict(pooled=np.nan, se=np.nan, tau2=np.nan, n=int(ok.sum()))
    m = S.meta(np.asarray(f)[ok], np.asarray(se)[ok])
    return dict(pooled=float(m["random"]), se=float(m["se_random"]), tau2=float(m["tau2"]), n=int(ok.sum()))


def stage_power(a):
    t0 = time.time()
    OUT.mkdir(exist_ok=True)
    st = setup(True)
    ent = st["ent"][:3] if a.check else st["ent"]
    nsets = 3 if a.check else len(st["sets"])
    U = st["U"]
    print(U[U.enters][["compartment", "family", "donors_400", "median_cells"]].to_string(index=False), flush=True)
    from multiprocessing import Pool as MP
    # (1) companions' psi_cells(400)
    jobs = [(k, j) for k in ent for j in range(nsets)]
    rows = []
    with MP(a.jobs, initializer=init_worker, initargs=(st,)) as mp:
        for i, r in enumerate(mp.imap_unordered(task_companion, jobs, chunksize=1)):
            rows.append(r)
            if (i + 1) % 50 == 0 or i + 1 == len(jobs):
                print(f"companions {i + 1}/{len(jobs)} {time.time() - t0:.0f}s", flush=True)
    C = pd.DataFrame(rows).sort_values(["compartment", "set"])
    C.to_csv(OUT / ("companions_c400_check.csv" if a.check else "companions_c400.csv"), index=False)
    lvl = C[C.valid >= N_DRAWS // 2].groupby("compartment").psi400.agg(P="median", q05=lambda x: x.quantile(.05),
                                                                         q95=lambda x: x.quantile(.95), n_sets="size")
    print(lvl.round(3).to_string(), flush=True)
    # (2) simulation
    reps, B = (1, 2) if a.check else (a.reps, a.boot)
    jobs = [(k, kf, rep, float(lvl.P[k]), B) for rep in range(reps) for kf in (0.0, 0.5, 1.0) for k in ent if k in lvl.index]
    sim = []
    last = time.time()
    with MP(a.jobs, initializer=init_worker, initargs=(st,)) as mp:
        for i, r in enumerate(mp.imap_unordered(task_sim, jobs, chunksize=1)):
            sim.append(r)
            if time.time() - last > 600 or i + 1 == len(jobs):
                last = time.time()
                print(f"sim {i + 1}/{len(jobs)} {time.time() - t0:.0f}s (last: {r['secs_point']:.0f}s point, "
                      f"{r['secs_boot']:.0f}s boot)", flush=True)
                pd.DataFrame(sim).to_csv(OUT / ("sim_check.csv" if a.check else "sim.csv"), index=False)
    Sim = pd.DataFrame(sim)
    Sim.to_csv(OUT / ("sim_check.csv" if a.check else "sim.csv"), index=False)
    rows = []
    for kf, g in Sim.groupby("kfac"):
        pooled = []
        for rep, h in g.groupby("rep"):
            m = pool_units(h.f, h.f_se)
            pooled.append(dict(rep=rep, pooled_f=m["pooled"], pooled_se=m["se"], units=m["n"],
                               pooled_f_alt=pool_units(h.f_alt, h.f_se)["pooled"]))
        pr = pd.DataFrame(pooled)
        ftrue = float(g.f_true.mean())
        rows.append(dict(kfac=kf, f_true=ftrue, units=g.compartment.nunique(), reps=pr.rep.nunique(),
                         unit_bias_mean=float((g.f - g.f_true).mean()), unit_bias_median=float((g.f - g.f_true).median()),
                         unit_f_sd=float(g.f.std()), unit_boot_se_median=float(g.f_se.median()),
                         pooled_f_mean=float(pr.pooled_f.mean()), pooled_bias=float(pr.pooled_f.mean() - ftrue),
                         pooled_se_boot_mean=float(pr.pooled_se.mean()), pooled_f_sd_over_reps=float(pr.pooled_f.std(ddof=1)),
                         pooled_bias_alt=float(pr.pooled_f_alt.mean() - ftrue), s2_over_true=float((g.s2 / g.s2_true.replace(0, np.nan)).median()),
                         psi400_median=float(g.psi400.median())))
    Pw = pd.DataFrame(rows)
    se_max = float(Pw.pooled_se_boot_mean.max())
    Pw["decidable"] = bool(ent and len(ent) >= 3 and se_max <= 0.15)
    Pw.to_csv(OUT / ("power_check.csv" if a.check else "power.csv"), index=False)
    print(Pw.round(4).T.to_string(), f"\nunits entering {len(ent)}, max pooled SE {se_max:.3f}, decidable "
          f"{bool(len(ent) >= 3 and se_max <= 0.15)}; runtime {time.time() - t0:.0f}s", flush=True)


# ---- score (implemented; run only after the power note) -----------------------------------------------------------

def task_score(job):
    """Clock or companion set j (-1 = clock) of one unit through the full pipeline; the clock also gets B donor bootstraps."""
    k, j, B = job
    pool, anc_k = G["pools"][k], anchors_without(k)
    if j < 0:
        X, genes, w, host = pool.clock, CLOCK, WF, host_rows_of(G["T"].pb, anc_k, CLOCK)
    else:
        w = G["sets"][j]
        cols = [G["cgenes"].index(g) for g in w.genes]
        X, genes, host = pool.comp_x[:, cols], list(w.genes), host_rows_of(G["gt"], anc_k, w.genes)
    ds, names = run_pipeline(pool, X, genes, host, w, anc_k)
    est = estimands(psi_of(ds, names))
    fb = boot_f(ds, names, B, np.random.default_rng([SEED, 11, pool.ci, j + 1])) if B else np.array([])
    return dict(compartment=k, set=j, **est, f_boot=";".join(f"{x:.5f}" for x in fb))


def stage_score(a):
    t0 = time.time()
    st = setup(True)
    from multiprocessing import Pool as MP
    jobs = [(k, j, a.boot if j < 0 else 0) for k in st["ent"] for j in range(-1, len(st["sets"]))]
    rows = []
    last = time.time()
    with MP(a.jobs, initializer=init_worker, initargs=(st,)) as mp:
        for i, r in enumerate(mp.imap_unordered(task_score, jobs, chunksize=1)):
            rows.append(r)
            if time.time() - last > 600 or i + 1 == len(jobs):
                last = time.time()
                print(f"score {i + 1}/{len(jobs)} {time.time() - t0:.0f}s", flush=True)
    R = pd.DataFrame(rows).sort_values(["compartment", "set"])
    R.to_csv(OUT / "score_units.csv", index=False)
    clk, cmp_ = R[R.set < 0].set_index("compartment"), R[R.set >= 0]
    fb = {k: np.array([float(x) for x in v.split(";")]) for k, v in clk.f_boot.items()}
    se = pd.Series({k: np.nanstd(v, ddof=1) for k, v in fb.items()})
    # H408.1: pooled clock s2 vs the q95 of pooled companion s2 (the clock bootstrap SEs weight both)
    s2se = pd.Series({k: np.nanstd(v * 400 * clk.psi400[k], ddof=1) for k, v in fb.items()})
    pooled_clock = pool_units(clk.s2, s2se)
    pooled_sets = np.array([pool_units(g.set_index("compartment").s2.reindex(clk.index), s2se)["pooled"]
                            for _, g in cmp_.groupby("set")])
    q95 = float(np.nanquantile(pooled_sets, 0.95))
    # H408.2: pooled f minus the companion-median f per unit, 90 % upper bound from the clock donor bootstraps
    cmed = cmp_.groupby("compartment").f.median().reindex(clk.index)
    diff = clk.f - cmed
    point = pool_units(diff, se)
    B = min(len(v) for v in fb.values())
    boots = np.array([pool_units(pd.Series({k: fb[k][b] for k in fb}) - cmed, se)["pooled"] for b in range(B)])
    upper = float(np.nanquantile(boots, 0.95))
    tests = dict(units=len(clk), s2_clock_pooled=pooled_clock["pooled"], s2_companion_q95=q95,
                 H408_1=bool(pooled_clock["pooled"] > q95), f_clock_minus_companion=point["pooled"],
                 f_upper90=upper, H408_2=bool(upper < 0.5), companion_f_median=float(cmed.median()))
    pd.DataFrame([tests]).to_csv(OUT / "tests.csv", index=False)
    print(pd.Series(tests).to_string(), f"\nruntime {time.time() - t0:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["prep", "power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--boot", type=int, default=40)
    a = ap.parse_args()
    limit_threads(1)
    warnings.filterwarnings("ignore")
    {"prep": stage_prep, "power": stage_power, "score": stage_score}[a.stage](a)


if __name__ == "__main__":
    main()
