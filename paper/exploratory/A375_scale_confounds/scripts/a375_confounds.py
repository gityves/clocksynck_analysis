"""A375 (registered): is the person's shared clock amplitude a batch, collection-centre or body-wide composition effect?

GTEx fast deaths, A336 rows (27 peripheral tissues). Per tissue pair, leverage = correlation of the two tissues' clock
radii after residualising on the phase harmonics plus a variant's covariates; L = median over pairs. Variants: none,
different isolation and expression batch only, collection centre, composition markers, genome PCs, all. Core loop
primary, 14 genes reported; 50 expression-matched companions.

    python a375_confounds.py
Outputs (`../output/`): tests.csv, pairs.csv. Aggregates only.
"""
import itertools
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A336_two_readouts" / "scripts"))
sys.path.insert(0, str(EXP / "A143_brain_artefact" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a336_person_gtex as P  # noqa: E402
import compute_brain_artefact_checks as A143  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE, RESPONSIVE_GENES  # noqa: E402
from clocksync._core.responsive_broad import RESPONSIVE_BROAD  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, pseudo_waveforms, states  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
MIN_DONORS, N_SETS = 40, 50
MARKERS = {"blood": ["HBB", "HBA1", "HBA2"], "immune": ["PTPRC", "LAPTM5", "CD53"], "endothelial": ["PECAM1", "VWF", "CDH5"],
           "fibroblast": ["COL1A1", "COL1A2", "DCN"], "smooth_muscle": ["ACTA2", "MYH11", "CNN1"]}
EXCL = set(CLOCK_LIKE) | set(RESPONSIVE_GENES) | set(RESPONSIVE_BROAD)
VARIANTS = ["none", "batch", "centre", "markers", "genome_pcs", "all"]


def say(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def attributes(pb):
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    o = a.obs[["subject", "tissue"]].copy()
    o["SAMPID"] = a.obs_names.astype(str)
    a.file.close()
    o["donor"], o["compartment"] = o.subject.astype(str), o.tissue.astype(str)
    o = o.drop_duplicates(["donor", "compartment"])
    at = pd.read_csv(A143.ATTRIBUTES, sep="\t", comment="#", low_memory=False,
                     usecols=["SAMPID", "SMNABTCH", "SMGEBTCH", "SMCENTER"])
    o = o.merge(at, on="SAMPID", how="left").set_index(["donor", "compartment"])
    return o.reindex(pd.MultiIndex.from_frame(pb[["donor", "compartment"]].astype(str))).reset_index(drop=True)


def tissue_covariates(x, depth, genes):
    """Per sample: five marker scores (z within tissue) and top 5 PCs of log CPM (2,000 most variable, exclusions)."""
    cpm = x / depth[:, None] * 1e6
    lg = np.log(cpm + 0.5)
    col = {g: k for k, g in enumerate(genes)}
    mk = []
    for gl in MARKERS.values():
        k = [col[g] for g in gl if g in col]
        s = lg[:, k].mean(1) if k else np.zeros(len(lg))
        mk.append((s - s.mean()) / (s.std() + 1e-12))
    keep = np.array([g not in EXCL for g in genes]) & (np.median(cpm, 0) >= 1)
    v = lg[:, keep]
    top = np.argsort(-v.var(0))[:2000]
    v = v[:, top] - v[:, top].mean(0)
    u, s, _ = np.linalg.svd(v, full_matrices=False)
    pcs = u[:, :5] * s[:5]
    pcs = (pcs - pcs.mean(0)) / (pcs.std(0) + 1e-12)
    return np.column_stack(mk), pcs


def lev(za, zb, ca, cb):
    ra, rb = np.abs(za), np.abs(zb)
    Hm = np.column_stack([np.ones(len(za))] + [g(k * p) for p in (np.angle(za), np.angle(zb)) for k in (1, 2) for g in (np.cos, np.sin)])
    res = lambda v, C: v - (X := np.column_stack([Hm, C]) if C is not None else Hm) @ np.linalg.lstsq(X, v, rcond=None)[0]  # noqa: E731
    xa, xb = res(ra, ca), res(rb, cb)
    return float(np.corrcoef(xa, xb)[0, 1])


def main():
    warnings.simplefilter("ignore")
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    pb, _ = P.load()
    pb["donor"] = pb.donor.astype(str)
    att = attributes(pb)
    genes = [c for c in pb.columns if c not in P.META]
    tissues = sorted(pb.compartment.unique())
    say(f"{len(tissues)} tissues, {pb.donor.nunique()} donors ({time.time() - t0:.0f}s)")
    x_all = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x_all / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x_all.mean(0) >= 1]
    lm = lm[~lm.index.duplicated()]
    T = {}
    for t in tissues:
        i = np.flatnonzero(pb.compartment.to_numpy() == t)
        mk, pcs = tissue_covariates(x_all[i], pb.depth.to_numpy(float)[i], genes)
        a = att.iloc[i]
        T[t] = dict(rows=i, donors=pb.donor.to_numpy()[i], mk=mk, pcs=pcs, na=a.SMNABTCH.astype(str).to_numpy(),
                    ge=a.SMGEBTCH.astype(str).to_numpy(), centre=(a.SMCENTER.astype(str).to_numpy() == "B1").astype(float))
    del x_all
    say("covariates done", f"({time.time() - t0:.0f}s)")
    sets = pseudo_waveforms(lm, N_SETS, np.random.default_rng([336, 2, 1]), P.WF, kind="matched")
    core_of = lambda w, j: Waveform(tuple(w.genes[k] for k in P.CORE_IDX), P.CORE.amp, P.CORE.psi, None, None, f"core{j}")  # noqa: E731
    wfs = [("core", "clock", P.CORE), ("14 genes", "clock", P.WF)] + [("core", f"matched{j}", core_of(w, j)) for j, w in enumerate(sets)]

    def zstates(w):
        out = {}
        for t in tissues:
            x = pb.iloc[T[t]["rows"]]
            L = projector(x[list(w.genes)], w)
            b, _ = states(x, L, w, x.depth.to_numpy(float))
            z = b[:, 0] + 1j * b[:, 1]
            out[t] = z - z.mean()
        return out

    pairs = list(itertools.combinations(tissues, 2))
    pos = {t: {d: k for k, d in enumerate(T[t]["donors"])} for t in tissues}
    rows = []
    for gs, lab, w in wfs:
        Z = zstates(w)
        for a, b in pairs:
            common = [d for d in T[a]["donors"] if d in pos[b]]
            if len(common) < MIN_DONORS:
                continue
            ia = np.array([pos[a][d] for d in common])
            ib = np.array([pos[b][d] for d in common])
            diff = (T[a]["na"][ia] != T[b]["na"][ib]) & (T[a]["ge"][ia] != T[b]["ge"][ib]) & (T[a]["na"][ia] != "nan")
            cen = T[a]["centre"][ia][:, None]
            spec = {"none": (None, None, slice(None)), "batch": (None, None, diff), "centre": (cen, cen, slice(None)),
                    "markers": (T[a]["mk"][ia], T[b]["mk"][ib], slice(None)),
                    "genome_pcs": (T[a]["pcs"][ia], T[b]["pcs"][ib], slice(None)),
                    "all": (np.column_stack([cen, T[a]["mk"][ia], T[a]["pcs"][ia]]), np.column_stack([cen, T[b]["mk"][ib], T[b]["pcs"][ib]]), diff)}
            for v, (ca, cb, m) in spec.items():
                za, zb = Z[a][ia][m], Z[b][ib][m]
                if len(za) < MIN_DONORS:
                    continue
                rows.append(dict(gene_set=gs, set=lab, variant=v, pair=f"{a}|{b}", donors=len(za),
                                 leverage=lev(za, zb, None if ca is None else ca[m], None if cb is None else cb[m])))
        if lab in ("clock",) or lab.endswith("9"):
            say(gs, lab, f"({time.time() - t0:.0f}s)")
    d = pd.DataFrame(rows)
    d[d.set == "clock"].to_csv(OUT / "pairs.csv", index=False)
    L = d.groupby(["gene_set", "set", "variant"]).leverage.median().reset_index()
    npairs = d.groupby(["gene_set", "set", "variant"]).pair.nunique().reset_index(name="pairs")
    L = L.merge(npairs)
    out = []
    for gs in ("core", "14 genes"):
        c = L[(L.gene_set == gs) & (L.set == "clock")].set_index("variant")
        for v in VARIANTS:
            r = dict(gene_set=gs, variant=v, L=float(c.loc[v, "leverage"]), pairs=int(c.loc[v, "pairs"]),
                     ratio_to_none=float(c.loc[v, "leverage"] / c.loc["none", "leverage"]))
            if gs == "core":
                m = L[(L.gene_set == "core") & (L.set.str.startswith("matched")) & (L.variant == v)].leverage
                r.update(matched_median=float(m.median()), matched_q95=float(m.quantile(.95)),
                         pct_matched=float((m < r["L"]).mean() * 100))
                if v != "none":
                    r["explains"] = bool(r["ratio_to_none"] < 0.5 or r["L"] <= r["matched_q95"])
            out.append(r)
    t = pd.DataFrame(out)
    t.to_csv(OUT / "tests.csv", index=False)
    print(t.round(4).to_string(index=False), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
