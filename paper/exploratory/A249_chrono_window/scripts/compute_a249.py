"""A249: are drug-relevant genes on body time in normal tissue and off it in the tumour of the same patient?

Registered in paper/docs/registered_experiments.md (A249). Body time is read from one sample of a patient and a gene from
the other. Blocks (a cancer of TCGA, or a single-cell cohort) hold, per patient: log2(CPM + 1) of the normal-side and of
the tumour-side sample (`yN`, `yT`) and the clock states of the two samples (`zN`, `zT`).
  S_N(set) = mean over genes of r(yN, zT)      S_T(set) = mean over genes of r(yT, zN)
r = |sum y conj(z)| / sqrt(sum y^2 sum |z|^2) per block (y, z centred within block or stratum), pooled over blocks with
patient-count weights (A229). Null: patients of y permuted within block, one permutation for all genes and both sides.
Companions: the state from 200 expression-matched pseudo-clock sets; 200 expression-matched gene sets as y.
D = S_N - S_T with patient bootstraps within block (states not refitted).

Designs
  tcga   the gated cancers of H138.4 (`q3_lib.tcga_states`), recount3 pairs; states on one projector per cancer (H138.4).
  sc     the A227 cohorts: normal epithelium vs the tumour-sample fibroblast state, tumour epithelium (all epithelial cells
         of the tumour sample; the genome-wide tables have no copy-number split) vs the normal-sample fibroblast state;
         states as `a149_lib.pair_states`.
  timed  GSE233242: mean over genes of the modulus of the two-term cosinor t against the recorded hour, normal minus tumour.
Stages
  power  eligibility, null SD and MDE (2.49 x SD) of S_N and S_T from permuted patients only; no observed association
                                                                 -> output/{power,genes}.csv
  score  (after the power note is dated) everything registered -> output/score_{tests,genes,matched,pseudo,descriptive}.csv
  `--placebo`: the patients of y are shuffled once within block before anything else (code test; output/placebo/).

    python compute_a249.py --stage power|score [--designs tcga sc timed] [--jobs 2] [--smoke]

Choices not fixed by the registration
- A gene's pooled r runs over the blocks where it is eligible (median count >= 10 in both sample types; TPM >= 1 in the
  timed cohort); a set is scored with >= half of its genes (>= 100 for F) eligible in >= 5 blocks, and S averages the
  genes eligible in >= 1 block.
- Pseudo-clock sets: drawn per block from the normal-side samples (TCGA: the normals; single cell: the normal-sample
  fibroblasts); set i is pooled across blocks. Matched gene sets: one gene per set gene from its 400 closest in mean
  normal-side log2(CPM + 1) (pooled over blocks), clock, waveform-like and set genes excluded.
- "Not scored" (power rule): < 100 pairs (TCGA), < 40 patients (single cell) or MDE > 0.10.
- Single cell: a block needs >= 4 patients with both fibroblast and both epithelial pseudobulks; y centred within stratum.
- Purity (descriptive): tertiles of ABSOLUTE purity within cancer; S recomputed in the low and the high tertile.
- Holm for H249.3: over sets K and F, separately for S_N and for S_T.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
for p in ("A227_tme_compartment_pairing", "A149_paired_tumour_normal", "A151_fibroblast_pairing", "A226_stroma_pairing_confirm",
          "A138_cancer_estimability"):
    sys.path.insert(0, str(PAPER / "exploratory" / p / "scripts"))
sys.path.insert(0, str(PAPER))
import clocksync as cs  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.paired import cohort_states, pseudo_waveforms  # noqa: E402
from common import cli  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parent / "output"
SEED = 249
WF = cs.gtex()
SETS = {"C": "DPYD TYMS TOP1 UGT1A1 CES1 CES2 ABCB1 ABCC2 CYP3A4 ERCC1 XPA BCL2 BAX TP53 WEE1".split(),
        "K": "WEE1 CCNB1 CDK1 MYC CCND1 GADD45A MDM2 CDKN1A MGMT".split(),
        "O": "BHLHE40 BHLHE41 NFIL3 HLF RORC ARNTL2".split()}
HYP = {"O": ("P249.0", "precondition"), "C": ("H249.1-2", "primary"), "K": ("H249.3", "secondary"), "F": ("H249.3", "secondary")}
N_PERM, N_PSEUDO, N_MATCH, N_POOL, N_BOOT, MDE_K = 10_000, 200, 200, 400, 2000, 2.49
MIN_COUNT, MIN_BLOCKS, MAX_MDE = 10, 5, 0.10
MIN_N = {"tcga": 100, "sc": 40}


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def fda_targets() -> list[str]:
    t = pd.read_csv(DATA / "drug_targets" / "raw" / "hpa_fda_targets.tsv", sep="\t")
    return sorted(set(t["Gene"].astype(str)))


def centre(v, stratum):
    return v - pd.DataFrame(v).groupby(stratum).transform("mean").to_numpy().reshape(np.shape(v))


def lcpm(counts: pd.DataFrame, depth) -> pd.DataFrame:
    return np.log2(counts.div(np.asarray(depth, float), axis=0) * 1e6 + 1).astype(np.float32)


def block(name, yN, yT, cN, cT, zN, zT, ZNp, ZTp, stratum=None):
    """Assemble a block: eligible genes, centred y matrices, centred states (and pseudo states)."""
    st = np.zeros(len(yN), int) if stratum is None else stratum
    ok = ((cN.median() >= MIN_COUNT) & (cT.median() >= MIN_COUNT)).to_numpy()
    genes = np.asarray(yN.columns)[ok]
    c = lambda z: centre(np.asarray(z), st)  # noqa: E731
    return dict(name=name, n=len(yN), genes=genes, stratum=st, mean_l=yN.loc[:, ok].mean().to_numpy(float),
                yN=centre(yN.loc[:, ok].to_numpy(float), st), yT=centre(yT.loc[:, ok].to_numpy(float), st),
                zN=c(zN), zT=c(zT), ZNp=np.array([c(z) for z in ZNp]), ZTp=np.array([c(z) for z in ZTp]))


# ---- designs ----------------------------------------------------------------------------------------------------------

def tcga_blocks(n_pseudo, smoke):
    import q3_lib as Q
    from common.datasets.tcga_paired import build as B
    gate = Q.tcga_states()
    out = []
    for i, c in enumerate([k for k, v in gate.items() if v["gated"]][:2 if smoke else None]):
        t, n, ct, cn = B._pairs(c)
        ct.index, cn.index = t.patient.to_numpy(), n.patient.to_numpy()
        ft, fn = ct.assign(depth=ct.sum(1)), cn.assign(depth=cn.sum(1))
        bt, _, bn, _ = cohort_states(ft, fn, True, WF)
        zT, zN = S.complex_state(bt), S.complex_state(bn)
        assert len(zT) == gate[c]["n"] and np.allclose(np.sort(abs(zT - zT.mean())), np.sort(abs(gate[c]["za"])), atol=1e-6), c
        frac = cn.div(cn.sum(1), axis=0)
        lm = np.log(frac.mean() + 1e-12)[(cn.mean() >= 1) | cn.columns.isin(list(WF.genes))]
        ws = pseudo_waveforms(lm, n_pseudo, np.random.default_rng([SEED, 1, i]), WF)
        ps = [cohort_states(ft, fn, True, w) for w in ws]
        out.append(block(c, lcpm(cn, fn.depth), lcpm(ct, ft.depth), cn, ct, zN, zT,
                         [S.complex_state(p[2]) for p in ps], [S.complex_state(p[0]) for p in ps]))
        say("tcga", c, len(zT), "pairs", len(out[-1]["genes"]), "genes")
    return out


def sc_blocks(n_pseudo, smoke):
    import a149_lib as L
    import compute_a227 as A
    fa, fb = "tumour:fibroblast", "adjacent_normal:fibroblast"
    out = []
    for i, s in enumerate(A.cohort_list(smoke)):
        pb = A.table(s)
        comps = set(pb.compartment)
        tum = L.MAL if L.MAL in comps else L.EPI
        if not {fa, fb, tum, L.NORM} <= comps:
            continue
        pts = sorted(set.intersection(*[set(pb[pb.compartment == c].donor) for c in (fa, fb, tum, L.NORM)]))
        u = L.pair_states(pb, fa, fb, pts=pts)
        if u is None:
            continue
        pts = u["pts"]
        G = A.genome_pair(s, (L.EPI, L.NORM, fa, fb), pts)   # genome tables: all tumour epithelium, no copy-number split
        genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth")]
        get = lambda c: G[G.compartment == c].drop_duplicates("donor").set_index("donor").reindex(pts)  # noqa: E731
        gt, gn, gfb = get(L.EPI), get(L.NORM), get(fb)
        if gt[genes[0]].isna().any() or gn[genes[0]].isna().any():
            continue
        frac = gfb[genes].div(gfb.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(gfb[genes].mean() >= 1) | frac.columns.isin(list(WF.genes))]
        ps = [p for w in pseudo_waveforms(lm, n_pseudo, np.random.default_rng([SEED, 2, i]), WF)
              if (p := L.pair_states(G, fa, fb, wf=w, pts=pts))]
        out.append(block(s, lcpm(gn[genes], gn.depth), lcpm(gt[genes], gt.depth), gn[genes], gt[genes], u["zb"], u["za"],
                         [p["zb"] for p in ps], [p["za"] for p in ps], stratum=np.asarray(u["stratum"])))
        say("sc", s, len(pts), "patients", len(out[-1]["genes"]), "genes")
    return out


# ---- statistics -------------------------------------------------------------------------------------------------------

def unit(Y):
    n = np.sqrt((Y ** 2).sum(0))
    return Y / np.where(n > 0, n, 1)


def r_obs(Y, z):
    """r of every column of Y (centred) with the state z (centred)."""
    return np.abs(unit(Y).T @ np.conj(z)) / np.sqrt(np.sum(abs(z) ** 2))


def index(blocks, genes):
    """Per block the column positions of `genes` (-1 where not eligible)."""
    return [pd.Index(b["genes"]).get_indexer(genes) for b in blocks]


def pooled(blocks, per_block, idx):
    """Pooled r per gene from per-block arrays (... x eligible genes of the block): weights = patients, over eligible blocks."""
    num = den = 0
    for b, r, ix in zip(blocks, per_block, idx):
        m = ix >= 0
        v = np.zeros(r.shape[:-1] + (len(ix),))
        v[..., m] = r[..., ix[m]]
        num, den = num + b["n"] * v, den + b["n"] * m
    return num / np.where(den > 0, den, np.nan)


def set_stat(blocks, genes, side, state_key=None, perms=None, rows=None):
    """S of a gene list on one side ('N': yN with zT; 'T': yT with zN). perms: per block n_perm x n permutations of y
    (returned as n_perm values); state_key: pseudo states ('ZTp' / 'ZNp') -> one S per pseudo set; rows: per block row
    index (bootstrap)."""
    idx = index(blocks, genes)
    per = []
    for k, (b, ix) in enumerate(zip(blocks, idx)):
        Y = b["y" + side][:, ix[ix >= 0]]
        z = b["zT" if side == "N" else "zN"]
        if rows is not None:
            Y, z = centre(Y[rows[k]], b["stratum"][rows[k]]), centre(z[rows[k]], b["stratum"][rows[k]])
        if state_key:
            Z = b[state_key]
            r = np.abs(np.conj(Z) @ unit(Y)) / np.sqrt((abs(Z) ** 2).sum(1))[:, None] if len(Z) else np.zeros((0, Y.shape[1]))
        elif perms is not None:
            r = np.abs(np.conj(z)[perms[k]] @ unit(Y)) / np.sqrt(np.sum(abs(z) ** 2))   # z permuted = y permuted
        else:
            r = r_obs(Y, z)
        out = np.zeros(r.shape[:-1] + (len(b["genes"]),))
        out[..., ix[ix >= 0]] = r
        per.append(out)
    if state_key:
        k = min(len(p) for p in per)
        per = [p[:k] for p in per]
    return np.nanmean(pooled(blocks, per, idx), axis=-1)


def eligible(blocks, genes):
    cnt = np.sum([ix >= 0 for ix in index(blocks, genes)], axis=0)
    return int((cnt >= MIN_BLOCKS).sum()), int((cnt >= 1).sum())


def matched_sets(blocks, genes, excl, rng):
    allg = sorted(set(np.concatenate([b["genes"] for b in blocks])) - set(excl))
    idx = index(blocks, allg)
    m = pooled(blocks, [b["mean_l"] for b in blocks], idx)
    ok = np.flatnonzero(np.sum([ix >= 0 for ix in idx], axis=0) >= min(MIN_BLOCKS, len(blocks)))
    tg = pooled(blocks, [b["mean_l"] for b in blocks], index(blocks, genes))
    picks = [rng.choice(ok[np.argsort(np.abs(m[ok] - t))[:N_POOL]], N_MATCH, replace=False) for t in tg if np.isfinite(t)]
    return [[allg[p[i]] for p in picks] for i in range(N_MATCH)]


def permutations(blocks, n_perm, seed):
    rng = np.random.default_rng(seed)
    out = []
    for b in blocks:
        P = np.tile(np.arange(b["n"]), (n_perm, 1))
        for s in np.unique(b["stratum"]):
            i = np.flatnonzero(b["stratum"] == s)
            P[:, i] = i[np.argsort(rng.random((n_perm, len(i))), axis=1)]
        out.append(P)
    return out


def chunked_null(blocks, genes, side, perms, size=1000):
    return np.concatenate([set_stat(blocks, genes, side, perms=[p[a:a + size] for p in perms])
                           for a in range(0, len(perms[0]), size)])


# ---- timed cohort -----------------------------------------------------------------------------------------------------

def cosinor_modulus(Y, tm):
    """|t_cos + i t_sin| of the two-term cosinor (intercept only) for every column of Y; tm: perms x n hours."""
    w = 2 * np.pi / 24
    Y = Y - Y.mean(0)
    c, s = np.cos(w * tm), np.sin(w * tm)
    c, s = c - c.mean(1, keepdims=True), s - s.mean(1, keepdims=True)
    cc, ss, sc = (c * c).sum(1)[:, None], (s * s).sum(1)[:, None], (c * s).sum(1)[:, None]
    det, cy, sy = cc * ss - sc ** 2, c @ Y, s @ Y
    bc, bs = (ss * cy - sc * sy) / det, (cc * sy - sc * cy) / det
    sig2 = ((Y ** 2).sum(0)[None, :] - bc * cy - bs * sy) / (len(Y) - 3)
    return np.hypot(bc / np.sqrt(sig2 * ss / det), bs / np.sqrt(sig2 * cc / det))


def timed(sets, n_perm, placebo):
    from common.datasets.breast_timed_gse233242 import build as B
    k = pd.read_csv(B.TRUTH)
    k = k[k.note.fillna("") == ""]
    x = pd.read_csv(B.TPM).drop_duplicates("gene_name").set_index("gene_name")
    h = k.resection_h.to_numpy(float)
    if placebo:
        h = h[np.random.default_rng([SEED, 9]).permutation(len(h))]
    P = np.argsort(np.random.default_rng([SEED, 5]).random((n_perm, len(h))), axis=1)
    rows = []
    for name, genes in sets.items():
        g = [a for a in genes if a in x.index]
        tn, tt = x.loc[g, k.tpm_column_normal].T, x.loc[g, k.tpm_column_tumour].T
        ok = ((tn.median() >= 1) & (tt.median() >= 1)).to_numpy()
        Y = {s: np.log2(t.loc[:, ok].to_numpy(float) + 1) for s, t in (("N", tn), ("T", tt))}
        mod = lambda s, tm: cosinor_modulus(Y[s], tm).mean(1)  # noqa: E731
        nul = mod("N", h[P]) - mod("T", h[P])
        rows.append(dict(set=name, genes=int(ok.sum()), patients=len(h), null_sd=float(nul.std()),
                         _obs=float(mod("N", h[None, :])[0] - mod("T", h[None, :])[0]),
                         _mN=float(mod("N", h[None, :])[0]), _mT=float(mod("T", h[None, :])[0]), _nul=nul))
    return rows


# ---- stages -----------------------------------------------------------------------------------------------------------

def run(args):
    sets = dict(SETS, F=fda_targets())
    excl = set(WF.genes) | {g for v in sets.values() for g in v}
    n_perm = 200 if args.smoke else N_PERM
    n_pseudo = 5 if args.smoke else N_PSEUDO
    power, tests, gene_rows, mrows, prows, desc = [], [], [], [], [], []
    for design in [d for d in args.designs if d != "timed"]:
        blocks = (tcga_blocks if design == "tcga" else sc_blocks)(n_pseudo if args.stage == "score" else 0, args.smoke)
        if args.placebo:
            rng = np.random.default_rng([SEED, 8])
            for b in blocks:
                p = rng.permutation(b["n"])
                b["yN"], b["yT"] = centre(b["yN"][p], b["stratum"]), centre(b["yT"][p], b["stratum"])
        perms = permutations(blocks, n_perm, [SEED, 3, len(blocks)])
        n_tot = sum(b["n"] for b in blocks)
        nulls = {}
        for name, genes in sets.items():
            e5, e1 = eligible(blocks, genes)
            need = 100 if name == "F" else len(genes) / 2
            row = dict(design=design, set=name, hypothesis=HYP[name][0], patients=n_tot, blocks=len(blocks), genes=len(genes),
                       genes_in_5_blocks=e5, genes_scored=e1)
            for side in "NT":
                nulls[(name, side)] = nul = chunked_null(blocks, genes, side, perms)
                row.update({f"null_mean_S_{side}": float(nul.mean()), f"null_sd_S_{side}": float(nul.std()), f"mde_S_{side}": MDE_K * float(nul.std())})
            row["scored"] = bool(e5 >= need and n_tot >= MIN_N[design] and max(row["mde_S_N"], row["mde_S_T"]) <= MAX_MDE)
            power.append(row)
        say(design, "nulls done")
        if args.stage == "power":
            continue
        for name, genes in sets.items():
            ms = matched_sets(blocks, genes, excl, np.random.default_rng([SEED, 4, list(sets).index(name)]))
            boot = []
            rng = np.random.default_rng([SEED, 7, list(sets).index(name)])
            for _ in range(50 if args.smoke else N_BOOT):
                rows = []
                for b in blocks:
                    r = np.arange(b["n"])
                    for s in np.unique(b["stratum"]):
                        i = np.flatnonzero(b["stratum"] == s)
                        r[i] = rng.choice(i, len(i))
                    rows.append(r)
                boot.append(set_stat(blocks, genes, "N", rows=rows) - set_stat(blocks, genes, "T", rows=rows))
            boot = np.array(boot)
            res = {}
            for side in "NT":
                s, nul = float(set_stat(blocks, genes, side)), nulls[(name, side)]
                pse = set_stat(blocks, genes, side, state_key="ZTp" if side == "N" else "ZNp")
                mat = np.array([set_stat(blocks, m, side) for m in ms])
                res[side] = dict(S=s, null_mean=float(nul.mean()), null_sd=float(nul.std()), p=float((1 + np.sum(nul >= s)) / (1 + len(nul))),
                                 pseudo_median=float(np.median(pse)), pseudo_q95=float(np.quantile(pse, .95)),
                                 p_clock=float((1 + np.sum(pse >= s)) / (1 + len(pse))), matched_median=float(np.median(mat)),
                                 matched_q95=float(np.quantile(mat, .95)), p_matched=float((1 + np.sum(mat >= s)) / (1 + len(mat))),
                                 mde=MDE_K * float(nul.std()))
                prows += [dict(design=design, set=name, side=side, i=i, value=float(v)) for i, v in enumerate(pse)]
                mrows += [dict(design=design, set=name, side=side, i=i, value=float(v)) for i, v in enumerate(mat)]
            pw = next(p for p in power if p["design"] == design and p["set"] == name)
            tests.append(dict(design=design, set=name, hypothesis=HYP[name][0], role=HYP[name][1], scored=pw["scored"], patients=n_tot,
                              blocks=len(blocks), genes_scored=pw["genes_scored"],
                              **{f"{k}_{sd}": v for sd in "NT" for k, v in res[sd].items()},
                              D=res["N"]["S"] - res["T"]["S"], D_lo90=float(np.quantile(boot, .05)), D_hi90=float(np.quantile(boot, .95)),
                              D_boot_median=float(np.median(boot))))
            idx = index(blocks, genes)
            for side in "NT":
                per = [r_obs(b["y" + side], b["zT" if side == "N" else "zN"]) for b in blocks]
                cr = [unit(b["y" + side]).T @ np.conj(b["zT" if side == "N" else "zN"]) / np.sqrt(np.sum(abs(b["zT" if side == "N" else "zN"]) ** 2))
                      for b in blocks]
                pr, pc = pooled(blocks, per, idx), pooled(blocks, [np.real(c) for c in cr], idx) + 1j * pooled(blocks, [np.imag(c) for c in cr], idx)
                nb = np.sum([ix >= 0 for ix in idx], axis=0)
                gene_rows += [dict(design=design, set=name, side=side, gene=g, blocks=int(k), r=float(r),
                                   implied_peak_h=float((-np.angle(c)) % (2 * np.pi) * S.HOURS) if k else np.nan)
                              for g, r, c, k in zip(genes, pr, pc, nb)]
            say(design, name, "scored")
        # descriptives: same-sample readouts, set O in the tumour is S_T above; blocks with >= 30 patients
        for name, genes in sets.items():
            idx = index(blocks, genes)
            for label, ykey, zkey in (("same_sample_normal", "yN", "zN"), ("same_sample_tumour", "yT", "zT")):
                per = [r_obs(b[ykey], b[zkey]) for b in blocks]
                desc.append(dict(design=design, item=label, set=name, value=float(np.nanmean(pooled(blocks, per, idx)))))
            for b in blocks:
                if b["n"] >= 30:
                    for side in "NT":
                        desc.append(dict(design=design, item=f"S_{side}", set=name, block=b["name"], patients=b["n"],
                                         value=float(set_stat([b], genes, side)), pseudo_median=float(np.median(
                                             set_stat([b], genes, side, state_key="ZTp" if side == "N" else "ZNp")))))
    if "timed" in args.designs:
        for r in timed(sets, n_perm, args.placebo):
            nul = r.pop("_nul")
            o, mN, mT = r.pop("_obs"), r.pop("_mN"), r.pop("_mT")
            power.append(dict(design="timed", hypothesis="H249.5", mde=MDE_K * r["null_sd"], **r))
            if args.stage == "score":
                desc.append(dict(design="timed", item="H249.5 normal minus tumour cosinor modulus", set=r["set"], value=o,
                                 modulus_normal=mN, modulus_tumour=mT, null_sd=r["null_sd"],
                                 p_one_sided=float((1 + np.sum(nul >= o)) / (1 + len(nul)))))
    pd.DataFrame(power).to_csv(OUT / "power.csv", index=False)
    print(pd.DataFrame(power).round(4).to_string(), flush=True)
    if args.stage == "score":
        T = pd.DataFrame(tests)
        T.to_csv(OUT / "score_tests.csv", index=False)
        pd.DataFrame(gene_rows).to_csv(OUT / "score_genes.csv", index=False)
        pd.DataFrame(mrows).to_csv(OUT / "score_matched.csv", index=False)
        pd.DataFrame(prows).to_csv(OUT / "score_pseudo.csv", index=False)
        pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
        print(T.round(4).T.to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--designs", nargs="+", default=["tcga", "sc", "timed"], choices=["tcga", "sc", "timed"])
    ap.add_argument("--placebo", action="store_true")
    ap.add_argument("--smoke", action="store_true", help="2 cancers / 3 cohorts, 200 permutations, 5 pseudo sets, 50 bootstraps")
    args = ap.parse_args(argv)
    cli.setup(args)
    global OUT
    OUT = OUT / ("placebo" if args.placebo else "smoke" if args.smoke else "")
    OUT.mkdir(parents=True, exist_ok=True)
    run(args)


if __name__ == "__main__":
    main()
