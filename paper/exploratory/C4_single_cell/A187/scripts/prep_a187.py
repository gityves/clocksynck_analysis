"""A187 prep: per-patient x bin (NC / S / G2M) pseudobulks of malignant cells from raw counts joined to the CoPhaser
cell-cycle tables. Reads counts with the loaders of the pilot fit script (the exec into the CoPhaser env is disabled; the
loaders only read). Run: python -I prep_a187.py [study ...]

Outputs (../output): pseudobulk_clock.csv.gz (study, patient, bin, n_cells, depth, 14 clock genes), genome_<study>.csv.gz
(genome-wide summed counts per patient x bin; pseudo-clock sets and pseudo_means), match.csv (join rates).
"""
import importlib.util, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp

HERE = Path(__file__).resolve()
C4 = HERE.parents[2]
OUT = HERE.parents[1] / "output"
STUDIES = ["GSE178341", "GSE132465", "GSE131907", "GSE176078", "GSE161529", "SRP327758", "GSE181919", "SRP416180"]
PILOT = {"GSE178341", "GSE132465", "GSE131907", "GSE176078", "SRP327758"}
MIN_UMI = 1000

_ex = os.execv
os.execv = lambda *a, **k: None  # the fit script switches to the CoPhaser env on import; loaders need no torch
spec = importlib.util.spec_from_file_location("fitcc", C4 / "pilot/scripts/fit_cellcycle.py")
F = importlib.util.module_from_spec(spec); sys.modules["fitcc"] = F
spec.loader.exec_module(F)
os.execv = _ex
from clocksync import paired as P
from clocksync.genes import canonical
CLOCK = list(P.gtex().genes)


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def cells_table(study):
    d = C4 / ("pilot/output" if study in PILOT else "fits/output")
    t = pd.read_csv(d / f"cells_{study}.csv.gz", index_col=0)
    assert t.index.is_unique
    return t


def load_counts(study):
    rng = np.random.default_rng(0)
    if study == "GSE176078":
        a = F.load_wu(3000, rng)
    elif study == "GSE161529":
        a = F.load_pal(3000, rng)
    elif (F.D / "scbc_full" / study).is_dir():
        a = F.load_srp(study, 3000, rng)
    else:
        a = F.load_geo(study, 3000, rng)
    return a


def process(study):
    t = cells_table(study)
    sp_, gp_ = (pd.read_csv(C4 / "fits/output/qc.csv").set_index("study").loc[study, ["S_peak_theta", "G2M_peak_theta"]])
    a = load_counts(study)
    X = sp.csr_matrix(a.X)
    names = pd.Index(a.obs_names)
    pos = names.get_indexer(t.index)
    ok = pos >= 0
    umi = np.asarray(X.sum(1)).ravel()
    umi_ok = np.isclose(umi[pos[ok]], t.total_umi.to_numpy()[ok])
    say(study, f"table cells {len(t)}, matched {ok.sum()} ({ok.mean():.4f}), total_umi equal {umi_ok.mean():.4f}")
    t = t[ok].copy(); X = X[pos[ok]]
    cyc = t.cycling_prob.to_numpy() >= 0.5
    d = lambda th, pk: np.abs(np.angle(np.exp(1j * (th - pk))))
    near_s = d(t.theta.to_numpy(), sp_) <= d(t.theta.to_numpy(), gp_)
    t["bin"] = np.where(~cyc, "NC", np.where(near_s, "S", "G2M"))
    gene = np.array([canonical(g) for g in a.var_names])
    cols = {g: np.flatnonzero(gene == g) for g in CLOCK}
    key = t.patient.astype(str) + "\t" + t["bin"]
    ukey, inv = np.unique(key.to_numpy(), return_inverse=True)
    M = sp.csr_matrix((np.ones(len(inv)), (inv, np.arange(len(inv)))), shape=(len(ukey), len(inv)))
    S = (M @ X).tocsr()
    rows = []
    for i, k in enumerate(ukey):
        p, b = k.split("\t")
        r = dict(study=study, patient=p, bin=b, n_cells=int(M[i].sum()), depth=float(t.total_umi.to_numpy()[inv == i].sum()))
        for g in CLOCK:
            r[g] = float(S[i][:, cols[g]].sum()) if len(cols[g]) else 0.0
        rows.append(r)
    pb = pd.DataFrame(rows)
    G = pd.DataFrame(S.toarray().astype(np.int32), index=pd.MultiIndex.from_tuples([tuple(k.split("\t")) for k in ukey], names=["patient", "bin"]),
                     columns=gene)
    G = G.T.groupby(level=0).sum().T
    G = G.loc[:, G.sum(0) > 0]
    G.to_csv(OUT / f"genome_{study}.csv.gz")
    m = dict(study=study, table_cells=len(ok), matched=int(ok.sum()), match_rate=float(ok.mean()), umi_equal_rate=float(umi_ok.mean()),
             malignant_cells_loaded=a.n_obs, S_peak=float(sp_), G2M_peak=float(gp_),
             **{f"cells_{b}": int((t["bin"] == b).sum()) for b in ("NC", "S", "G2M")})
    return pb, m


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    todo = sys.argv[1:] or STUDIES
    pbs, ms = [], []
    for s in todo:
        pb, m = process(s)
        pb.to_csv(OUT / f"pseudobulk_clock_{s}.csv.gz", index=False)
        pd.DataFrame([m]).to_csv(OUT / f"match_{s}.csv", index=False)
        del pb
    say("done")
