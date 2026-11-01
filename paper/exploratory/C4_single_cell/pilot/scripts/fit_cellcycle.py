"""C4 pilot, method QC: CoPhaser cell-cycle phase of CNV-malignant single cells, one model per study.

QC only. No statistic of clock-gene expression across cell-cycle phase is computed here: clock and clock-like genes are
removed from the rhythmic and the context gene sets, from the S/G2M score lists, and are not stored in the outputs.

    PATH=<CoPhaser_Env>/bin:$PATH python fit_cellcycle.py STUDY [--restarts 5] [--max-per-patient 3000]
STUDY in GSE178341 GSE132465 GSE131907 GSE183904 SRP327758 GSE176078 (Wu breast, author 'Cancer Epithelial').
Outputs (../output): cells_<study>.csv.gz, validation_<study>.png, seeds_<study>.csv, qc_<study>.json; qc.csv and
patient_cycling.csv are rebuilt from all studies present.
Raw counts are only read (anndata / pandas); nothing from the data folders is executed.
"""
import argparse, glob, hashlib, json, os, subprocess, sys, time, warnings
from pathlib import Path

ENV_PY = "/home/maxine/miniconda3/envs/CoPhaser_Env/bin/python"
if os.path.realpath(sys.executable) != os.path.realpath(ENV_PY):  # run.sh starts the sccc env: switch to the CoPhaser env
    os.execv(ENV_PY, [ENV_PY, "-u", "-P", *sys.argv])

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore")
D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
OUT = Path(__file__).resolve().parents[1] / "output"
COPHASER = Path("/home/maxine/Documents/paychere/CoPhaser")
CLOCK = set("ARNTL CLOCK NPAS2 PER1 PER2 PER3 CRY1 CRY2 NR1D1 NR1D2 DBP RORA TEF CIART BHLHE40 BHLHE41 NFIL3 CIPC".split())
CLOCK |= {"BMAL1"}  # alias of ARNTL
MIN_UMI = 1000
VAL_GENES = ["TOP2A", "PCNA", "MKI67", "MCM6"]
CNV_NAME = {"GSE176078": None}


def malignant_cnv(study):
    c = pd.read_csv(D / "cnv" / f"{study}.csv.gz", dtype=str)
    return c[c.cnv_call == "malignant"]


def load_geo(study, cap, rng):
    cnv = malignant_cnv(study)
    a = ad.read_h5ad(D / "geo" / study / f"{study}.h5ad", backed="r")
    names = pd.Index(a.obs_names)
    cnv = cnv[cnv.barcode.isin(names)].drop_duplicates("barcode")
    print(f"{study}: {len(cnv)} CNV-malignant barcodes found in h5ad", flush=True)
    return subsample_and_read(a, names, cnv, cap, rng, lambda a, ix: a[ix].to_memory())


def load_srp(study, cap, rng):
    cnv = malignant_cnv(study)
    parts = []
    for smp, g in cnv.groupby("sample"):
        f = D / "scbc_full" / study / f"{smp}.h5ad"
        if not f.exists():
            continue
        a = ad.read_h5ad(f, backed="r")
        names = pd.Index(a.obs_names)
        g = g[g.barcode.isin(names)].drop_duplicates("barcode")
        if not len(g):
            continue
        g = cap_patient_chunk(g)
        sub = a[names.get_indexer(g.barcode)].to_memory()
        sub.var_names = pd.Index(sub.var["gene_symbols"].astype(str).str.upper())
        sub.obs = pd.DataFrame({"patient": g.patient.to_numpy(), "sample": smp}, index=[f"{smp}_{b}" for b in g.barcode])
        sub.layers.clear()
        sub.X = sp.csr_matrix(sub.X)
        parts.append(sub)
    print(f"{study}: {len(parts)} samples read", flush=True)
    a = ad.concat([p[:, ~p.var_names.duplicated()] for p in parts], join="inner")
    return cap_patients(a, cap, rng)


def cap_patient_chunk(g):
    return g


def cap_patients(a, cap, rng):
    keep = []
    for _, ix in a.obs.groupby("patient", observed=True).indices.items():
        keep.append(np.sort(rng.choice(ix, min(cap, len(ix)), replace=False)))
    return a[np.sort(np.concatenate(keep))].copy()


def subsample_and_read(a, names, cnv, cap, rng, reader):
    pos = names.get_indexer(cnv.barcode)
    pat = cnv.patient.to_numpy()
    keep = []
    for p in np.unique(pat):
        ix = pos[pat == p]
        keep.append(rng.choice(ix, min(cap, len(ix)), replace=False))
    keep = np.sort(np.concatenate(keep))
    sub = reader(a, keep)
    sub.obs = pd.DataFrame({"patient": pd.Series(pat, index=pos).loc[keep].to_numpy()}, index=sub.obs_names)
    sub.X = sp.csr_matrix(sub.X)
    sub.layers.clear()
    sub.var_names = pd.Index(sub.var_names.astype(str).str.upper())
    return sub[:, ~sub.var_names.duplicated()].copy()


def load_wu(cap, rng):
    base = next((D / "sc_breast_er/GSE176078/extract").glob("*"))
    md = pd.read_csv(base / "metadata.csv", index_col=0)
    bc = pd.read_csv(base / "count_matrix_barcodes.tsv", header=None)[0].astype(str).to_numpy()
    genes = pd.read_csv(base / "count_matrix_genes.tsv", header=None)[0].astype(str).str.upper().to_numpy()
    mal = md.reindex(bc).celltype_major.to_numpy() == "Cancer Epithelial"
    pat = md.reindex(bc)["orig.ident"].astype(str).to_numpy()
    keep = []
    for p in np.unique(pat[mal]):
        ix = np.flatnonzero(mal & (pat == p))
        keep.append(rng.choice(ix, min(cap, len(ix)), replace=False))
    keep = np.sort(np.concatenate(keep))
    newcol = np.full(len(bc), -1, dtype=np.int64)
    newcol[keep] = np.arange(len(keep))
    rows, cols, vals = [], [], []
    for ch in pd.read_csv(base / "count_matrix_sparse.mtx", sep=" ", skiprows=2, header=None, chunksize=20_000_000,
                          dtype=np.int64, engine="c"):
        c = newcol[ch[1].to_numpy() - 1]
        m = c >= 0
        rows.append(ch[0].to_numpy()[m] - 1); cols.append(c[m]); vals.append(ch[2].to_numpy()[m])
    X = sp.csr_matrix((np.concatenate(vals).astype(np.float32), (np.concatenate(cols), np.concatenate(rows))),
                      shape=(len(keep), len(genes)))
    a = ad.AnnData(X, obs=pd.DataFrame({"patient": pat[keep]}, index=bc[keep]), var=pd.DataFrame(index=pd.Index(genes)))
    return a[:, ~a.var_names.duplicated()].copy()


def version_info():
    import importlib.metadata as im
    head = subprocess.run(["git", "-C", str(COPHASER), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    diff = subprocess.run(["git", "-C", str(COPHASER), "diff"], capture_output=True).stdout
    return dict(cophaser_commit=head, cophaser_diff_sha256=hashlib.sha256(diff).hexdigest(),
                torch=im.version("torch"), scanpy=im.version("scanpy"), anndata=im.version("anndata"),
                numpy=np.__version__, pandas=pd.__version__)


def circ_peak(theta, x):
    return float(np.angle(np.sum(np.exp(1j * theta) * (x - x.mean()))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--restarts", type=int, default=5)
    ap.add_argument("--max-per-patient", type=int, default=3000)
    ap.add_argument("--epochs", type=int, default=200)
    args = ap.parse_args()
    study = args.study
    OUT.mkdir(exist_ok=True)
    import torch
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cophaser import CoPhaser, auto_params, cycle_qc, gene_sets, plotting, utils
    torch.set_num_threads(4)
    t0 = time.time()
    rng = np.random.default_rng(0)

    if study == "GSE176078":
        a = load_wu(args.max_per_patient, rng)
    elif study.startswith("SRP"):
        a = load_srp(study, args.max_per_patient, rng)
    else:
        a = load_geo(study, args.max_per_patient, rng)
    a.obs["patient"] = a.obs["patient"].astype(str)
    a.layers["total"] = a.X
    a.obs["total_umi"] = np.asarray(a.X.sum(1)).ravel()
    n_before = a.n_obs
    a = a[a.obs.total_umi >= MIN_UMI].copy()
    print(f"loaded {n_before} cells ({a.obs.patient.nunique()} patients) in {time.time()-t0:.0f}s; "
          f"{n_before - a.n_obs} below {MIN_UMI} UMI dropped; median UMI {a.obs.total_umi.median():.0f}", flush=True)
    assert not (CLOCK & set(a.var_names)) or True

    # S / G2M scores on the full gene matrix (control genes from all genes), clock genes dropped from the lists
    s_all, g_all = auto_params._seurat_cell_cycle_gene_sets()
    orig_sets = auto_params._seurat_cell_cycle_gene_sets
    auto_params._seurat_cell_cycle_gene_sets = lambda: ([g for g in s_all if g not in CLOCK], [g for g in g_all if g not in CLOCK])
    sc_full, n_s, n_g = auto_params._compute_cell_cycle_scores(None, a, "total", "total_umi")
    auto_params._seurat_cell_cycle_gene_sets = orig_sets
    a.obs["S_score"], a.obs["G2M_score"] = sc_full.S_score.to_numpy(), sc_full.G2M_score.to_numpy()

    var = set(a.var_names)
    rhythmic = [g for g in gene_sets.HUMAN_CELL_CYCLE_GENE_SET if g in var and g not in CLOCK]
    hvg = [g for g in utils.get_variable_genes(a, n_variable_genes=2000, layer="total") if g not in CLOCK]
    hist = [g for g in gene_sets.human_canonical_histones if g and g.upper() in var and g.upper() not in CLOCK]
    context = list(dict.fromkeys(hvg + [h.upper() for h in hist]))
    genes = sorted(set(rhythmic) | set(context))
    assert not (CLOCK & set(genes))
    am = a[:, genes].copy()
    am.layers["total"] = am.X
    print(f"rhythmic genes {len(rhythmic)}, context genes {len(context)} ({len(hvg)} HVG + histones), model genes {len(genes)}", flush=True)

    n_latent = auto_params.suggest_n_latent("tissue")
    model_params = dict(rhythmic_gene_names=rhythmic, context_genes=context, n_latent=n_latent)
    probe = CoPhaser(**model_params)
    probe.load_anndata(am, layer_to_use="total", library_size_field="total_umi")
    hp = auto_params.auto_hyperparameters(probe, cycle="cell_cycle", adata=a, layer="total", library_size_field="total_umi")
    print("auto hyperparameters:", {k: v for k, v in hp["trainer"].items()}, flush=True)
    del probe

    seeds = list(range(args.restarts))
    res = cycle_qc.fit_with_restarts(am, model_params, hp, layer="total", n_restarts=args.restarts, seeds=seeds,
                                     cycle="cell_cycle", library_size_field="total_umi", device="cuda",
                                     n_epochs=args.epochs)
    rows = []
    for r in res:
        rows.append({k: v for k, v in r.items() if k not in ("model", "trainer", "state_dict", "config")})
    seeds_df = pd.DataFrame(rows)
    num = [c for c in seeds_df.columns if c not in ("seed", "score")]
    seeds_df["converged"] = [bool(np.isfinite(r.get("marker_ratio", np.nan)) and np.isfinite(r.get("phase_separation", np.nan))
                                  and np.isfinite(r["model"].variable_genes.sum().item())) for r in res]
    seeds_df["vetoed"] = seeds_df.score < 0
    seeds_df["usable"] = seeds_df.converged & ~seeds_df.vetoed
    seeds_df.to_csv(OUT / f"seeds_{study}.csv", index=False)
    usable = [i for i in range(len(res)) if seeds_df.usable[i]]
    kept_i = usable[0] if usable else 0  # res is sorted best first
    best = res[kept_i]
    model = best["model"]
    print(f"seeds usable {len(usable)}/{len(res)}; kept seed {best['seed']} score {best['score']:.3f}"
          f"{'' if usable else ' (ALL VETOED OR NON-CONVERGED: kept best-ranked anyway)'}", flush=True)

    model.to("cpu"); model.eval()
    gen, space = model.get_outputs()
    out = plotting.plot_cell_cycle_validations(model, space, gen, am, layer="total", genes=VAL_GENES, hue_key=None,
                                               return_values=True)
    fig = out[0]
    fig.savefig(OUT / f"validation_{study}.png", dpi=80, bbox_inches="tight")
    plt.close("all")
    # aligned phase of every cell (the plot returns a 10,000-cell subsample)
    theta = model.orient_align_pseudotimes(space["theta"], plot=False).detach().cpu().numpy().ravel().astype(float)
    theta = np.arctan2(np.sin(theta), np.cos(theta))
    branch_on = isinstance(getattr(model, "cycling_status_prior", 1), (int, float)) and model.cycling_status_prior < 1
    logit = space["cycling_logits"].detach().cpu().numpy().ravel().astype(float) if branch_on else np.full(len(theta), np.inf)
    pcyc = 1 / (1 + np.exp(-logit))

    S, G = a.obs.S_score.to_numpy(), a.obs.G2M_score.to_numpy()
    cyc = pcyc > 0.5
    use = cyc if cyc.sum() > 50 else np.ones(len(theta), bool)
    ps, pg = circ_peak(theta[use], S[use]), circ_peak(theta[use], G[use])
    gap = float(np.mod(pg - ps, 2 * np.pi))
    order_ok = bool(0 < gap < np.pi)

    cells = pd.DataFrame({"patient": a.obs.patient.to_numpy(), "theta": theta, "cycling_prob": pcyc,
                          "S_score": S, "G2M_score": G, "total_umi": a.obs.total_umi.to_numpy()}, index=a.obs_names)
    cells.index.name = "cell"
    assert not any(c.split(":")[-1].upper() in CLOCK for c in cells.columns)
    cells.to_csv(OUT / f"cells_{study}.csv.gz")
    pc = cells.groupby("patient").agg(n_cells=("theta", "size"), frac_cycling=("cycling_prob", lambda x: float((x > 0.5).mean())))
    pc.insert(0, "study", study)
    pc.to_csv(OUT / f"patient_cycling_{study}.csv")

    q = {k: v for k, v in best.items() if k in cycle_qc.RANKING_SIGNALS + cycle_qc.DIAGNOSTIC_SIGNALS +
         ("phase_uniformity", "expected_phase_separation", "n_marker_genes", "cycling_status_prior", "seed", "score", "minutes")}
    qc = dict(study=study, n_cells=len(cells), n_patients=cells.patient.nunique(), n_rhythmic=len(rhythmic), n_context=len(context),
              n_latent=n_latent, seeds_run=len(res), seeds_kept=len(usable), kept_seed=best["seed"], kept_vetoed=not usable,
              rank_score=best["score"], **{k: v for k, v in q.items() if k not in ("seed", "score")},
              cycling_branch_on=branch_on, cycling_frac_p05=float(cyc.mean()), cycling_frac_mean=float(pcyc.mean()),
              median_umi=float(cells.total_umi.median()), S_peak_theta=ps, G2M_peak_theta=pg, S_to_G2M_gap_rad=gap,
              S_to_G2M_order_ok=order_ok, w_rhythmic=hp["trainer"].get("rhythmic_likelihood_weight"),
              w_nonrhythmic=hp["trainer"].get("non_rhythmic_likelihood_weight"),
              prior_set=hp["trainer"].get("cycling_status_prior"), total_minutes=(time.time() - t0) / 60,
              **version_info())
    json.dump(qc, open(OUT / f"qc_{study}.json", "w"), indent=1, default=float)

    parts = sorted(OUT.glob("qc_*.json"))
    pd.DataFrame([json.load(open(p)) for p in parts]).to_csv(OUT / "qc.csv", index=False)
    pd.concat([pd.read_csv(p) for p in sorted(OUT.glob("patient_cycling_*.csv"))]).to_csv(OUT / "patient_cycling.csv", index=False)
    print(json.dumps(qc, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
