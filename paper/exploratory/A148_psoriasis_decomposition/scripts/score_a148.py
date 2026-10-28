"""A148 scoring (registered 2026-10-28): keratinocyte state-mix decomposition of the psoriatic delay, healthy -> uninvolved -> lesional.

Steps as registered: (1) per-sample NNLS shares on CPM against the five keratinocyte bin references plus fibroblast, endothelial,
immune, melanocyte (genes: 500 most bin-variable, no clock / cue / locked genes), shares renormalised within keratinocytes;
ACCEPTANCE CHECK (lesional - uninvolved shares: cycling up, most-differentiated kc_d4 down, pooled) before any prediction;
(2)-(3) composition-only predicted step offset arg(sum g_k^A z_k / sum g_k^B z_k) (A69; g = A60 gain of the group's mean-share
mixture, z = A63 healthy per-bin complex state), per series and pooled with the A54 weights (arg of the weighted sum of unit
complex ratios); (4) interval: 2,000 resamples of samples within group and series, each paired with a draw (recycled) of the 500
A63 donor-bootstrap fits. H148.1: lesional-uninvolved < 0, H148.2: uninvolved-healthy > 0, Holm over the two; D148: observed
(A115 within-group core-minus-cue estimator, same samples) vs predicted. Seed 148.

    python score_a148.py [--jobs 4]          (needs output/ref_counts.csv.gz from build_refs.py)
Outputs (output/): acceptance.csv, shares.csv, h148.csv, d148.csv, z_boot.csv
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import nnls
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER.parent)); sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A115_A116_psoriasis/scripts"))
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("psoriasis"))); sys.path.insert(0, str(group_dir("psoriasis_details")))
import _lib as A  # noqa: E402  (A115)
import _attribution as M  # noqa: E402
import _singlecell as C  # noqa: E402
import _shared as H  # noqa: E402
from compute_differentiation import healthy_table  # noqa: E402
from compute_bin_mix import weights as bin_weights, zs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[0] / "output"
SEED, NB, NZ = 148, 2000, 500
HOURS = A.HOURS
KC = ["kc_cycling", "kc_d1", "kc_d2", "kc_d3", "kc_d4"]
NONKC = ["fibroblast", "endothelial", "immune", "melanocyte"]
SERIES = ["GSE13355", "GSE121212", "GSE201827", "GSE268958"]
GROUPS = ("healthy", "nonlesional", "lesional")
LESIONAL_SRC = {"GSE13355": "processed_heldout", "GSE121212": "processed", "GSE201827": "processed", "GSE268958": "processed"}
NSELECT = 500
META = set(A.META) | {"dataset", "n_replicates"}


def read_all(path):
    d = pd.read_parquet(path)
    return d.rename(columns={"BMAL1": "ARNTL"}) if "ARNTL" not in d and "BMAL1" in d else d


def load_series():
    """series -> {group: DataFrame (one row per patient, all genes)}; lesional from the original processed matrices."""
    out = {}
    for s in SERIES:
        a = read_all(A.ROOT / "processed_a115" / f"{s}_psoriasis.parquet")
        l = read_all(A.ROOT / LESIONAL_SRC[s] / f"{s}.parquet")
        l = l[(l.disease == "psoriasis") & (l.status == "lesional")]
        out[s] = {"healthy": a[a.status == "healthy"].drop_duplicates("patient").reset_index(drop=True),
                  "nonlesional": a[a.status == "nonlesional"].drop_duplicates("patient").reset_index(drop=True),
                  "lesional": l.drop_duplicates("patient").reset_index(drop=True)}
    return out


def cpm(d, vt, gene_mean_log2=None):
    """CPM (samples x genes). counts: count / depth * 1e6; arrays: 2^x renormalised to 1e6 over the genes present."""
    g = [c for c in d.columns if c not in META]
    x = d[g].astype(float)
    if vt == "counts":
        return x.fillna(0).div(d.depth.astype(float).to_numpy(), axis=0) * 1e6
    if vt == "log2_centred":                         # gene-wise centred as provided: levels imputed (sensitivity only)
        x = x.add(gene_mean_log2.reindex(g), axis=1)
    lin = 2.0 ** x
    return lin.div(lin.sum(1), axis=0) * 1e6


def pick_genes(ref, avail, excl, n=NSELECT, min_cpm=20.0, cols=KC):
    c = ref[cols].div(ref[cols].sum(0), axis=1) * 1e6
    c = c.loc[c.index.intersection(avail).difference(excl)]
    c = c[c.mean(1) >= min_cpm]
    v = np.log2(c + 1).var(1)
    return list(v.nlargest(n).index)


def shares_of(y, R):
    """NNLS coefficients (samples x refs) of CPM rows y (samples x genes) on R (genes x refs)."""
    A_ = R.to_numpy()
    co = np.array([nnls(A_, v)[0] for v in y.to_numpy()])
    return pd.DataFrame(co, index=y.index, columns=R.columns)


def kc_shares(co):
    k = co[KC]; t = k.sum(1)
    return k.div(t.where(t > 0), axis=0), co[KC].sum(1) / co.sum(1).where(co.sum(1) > 0)


def gain_fn(Lp, cpmk):
    """group mean KC shares (Series over KC) -> complex gain per KC bin (A60)."""
    return lambda m: M.gains(M.shares(m, cpmk), Lp)


def pred_offset(gA, gB, z):
    ks = KC
    return float(np.angle(sum(gA[k] * z[k] for k in ks) / sum(gB[k] * z[k] for k in ks)) * HOURS)


def pool_unit(ratios, w):
    """arg of the weighted sum of unit complex ratios (h)."""
    u = np.array([r / abs(r) for r in ratios])
    return float(np.angle(np.sum(np.asarray(w) * u)) * HOURS)


def ratio(gA, gB, z):
    return sum(gA[k] * z[k] for k in KC) / sum(gB[k] * z[k] for k in KC)


def holm(p):
    o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for r, i in enumerate(o):
        run = max(run, (m - r) * p[i]); adj[i] = min(1, run)
    return adj


def observed(items, sets, w54, keys, tag):
    """A115 within-group estimator: per series, group offset arg E[c_core c_cue*] with one projector per series over the two
    groups of the step; step = wrapped difference; pooled with A54 weights (arg of weighted sum of complex group statistics)."""
    res = {}
    for step, (a, b) in {"L-U": ("lesional", "nonlesional"), "U-H": ("nonlesional", "healthy")}.items():
        cA, cB = {}, {}
        for s in keys:
            x = pd.concat([items[s][a], items[s][b]], ignore_index=True)
            counts = bool((x.value_type == "counts").all()); ia = np.arange(len(items[s][a])); ib = len(items[s][a]) + np.arange(len(items[s][b]))
            bc, vc = A.states_for(x, A.WF, counts)[0]; bk, vk = A.states_for(x, sets["cue"], counts)[0]
            cA[s] = A.S.ccorr(bc[ia], bk[ia], vc[ia], vk[ia]); cB[s] = A.S.ccorr(bc[ib], bk[ib], vc[ib], vk[ib])
            res[(step, s)] = A.wrap((np.angle(cA[s]) - np.angle(cB[s])) * HOURS)
        w = np.array([w54.get(s, 1.0) for s in keys])
        res[(step, "pooled")] = A.wrap((np.angle(np.sum(w * np.array([cA[s] for s in keys]))) - np.angle(np.sum(w * np.array([cB[s] for s in keys])))) * HOURS)
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--variant", choices=("registered", "allrefs"), default="registered", help="allrefs: POST HOC gene selection over all nine references, acceptance only")
    args = ap.parse_args()
    limit_threads(4)
    rng = np.random.default_rng(SEED)
    sets = A.cue_sets(); w54 = A.weights54()
    ref = pd.read_csv(OUT / "ref_counts.csv.gz", index_col=0)
    refc = ref.div(ref.sum(0), axis=1) * 1e6                  # CPM over all genes of the reference
    items = load_series()
    VT = {s: items[s]["nonlesional"].value_type.iloc[0] for s in SERIES}
    meta = META
    avail = set.intersection(*[set(c for c in items[s]["nonlesional"].columns if c not in meta) for s in SERIES]) & set(ref.index)
    excl = set(A.CORE) | set(sets["cue"].genes) | set(sets["locked"].genes) | {"BMAL1", "C1orf51", "CHRONO"}
    post = args.variant != "registered"; sfx = "_posthoc" if post else ""
    genes = pick_genes(ref, avail, excl, cols=KC + NONKC if post else KC)
    pd.Series(genes, name="gene").to_csv(OUT / f"genes_selected{sfx}.csv", index=False)
    R = refc.loc[genes, KC + NONKC]
    print("genes", len(genes), flush=True)
    # level imputation for the gene-wise centred series (sensitivity): mean log2 of GSE201827 over all its samples
    g201 = items["GSE201827"]
    gm = pd.concat([g201[g] for g in GROUPS])[[c for c in g201["healthy"].columns if c not in meta]].astype(float).mean()
    # ---- shares
    rows, SH = [], {}
    for s in SERIES:
        for g in GROUPS:
            d = items[s][g]; y = cpm(d, VT[s], gm)[genes]
            co = shares_of(y, R); ks, kf = kc_shares(co)
            SH[(s, g)] = ks
            res = ((y.to_numpy() - co.to_numpy() @ R.to_numpy().T) ** 2).sum() / (y.to_numpy() ** 2).sum()
            for b in KC:
                rows.append(dict(series=s, group=g, bin=b, mean_share=ks[b].mean(), sd=ks[b].std(ddof=1), n=len(d)))
            for b in NONKC:
                rows.append(dict(series=s, group=g, bin=b, mean_share=(co[b] / co.sum(1)).mean(), sd=(co[b] / co.sum(1)).std(ddof=1), n=len(d), note="share of all references"))
            rows.append(dict(series=s, group=g, bin="kc_fraction_of_all", mean_share=kf.mean(), sd=kf.std(ddof=1), n=len(d)))
            rows.append(dict(series=s, group=g, bin="unexplained_fraction_of_ss", mean_share=res, sd=np.nan, n=len(d)))
    shares_tab = pd.DataFrame(rows)
    shares_tab.to_csv(OUT / f"shares{sfx}.csv", index=False)
    # ---- acceptance check
    def pooled_diff(keys, b, rr=None):
        w = np.array([w54.get(s, 1.0) for s in keys]); w = w / w.sum()
        out = []
        for s in keys:
            L_, U_ = SH[(s, "lesional")][b], SH[(s, "nonlesional")][b]
            if rr is not None:
                L_, U_ = rr.choice(L_.to_numpy(), len(L_)), rr.choice(U_.to_numpy(), len(U_))
            out.append(np.nanmean(L_) - np.nanmean(U_))
        return float(np.dot(w, out)), out
    acc = []
    for tagset, keys in (("primary_3series", [s for s in SERIES if s != "GSE13355"]), ("sensitivity_4series", SERIES)):
        for b in KC:
            d, per = pooled_diff(keys, b)
            bs = np.array([pooled_diff(keys, b, rng)[0] for _ in range(1000)])
            acc.append(dict(set=tagset, bin=b, pooled_L_minus_U=d, ci_lo=np.quantile(bs, .025), ci_hi=np.quantile(bs, .975),
                            **{f"{s}": v for s, v in zip(keys, per)}))
    acc = pd.DataFrame(acc)
    ok = {}
    for tagset in acc.set.unique():
        a = acc[acc.set == tagset].set_index("bin").pooled_L_minus_U
        ok[tagset] = bool(a["kc_cycling"] > 0 and a["kc_d4"] < 0)
    acc["acceptance_pass"] = acc.set.map(ok)
    acc.to_csv(OUT / f"acceptance{sfx}.csv", index=False)
    print(acc.round(4).to_string(), "\nACCEPTANCE:", ok, flush=True)
    if post:
        return
    if not ok["primary_3series"]:
        print("ACCEPTANCE FAILED (primary): stop, no predicted offset computed", flush=True); return
    # ---- A63 healthy states and donor bootstrap
    pbh, _ = C.load("healthy_skin_bins")
    t = healthy_table(pbh, ["healthy"])
    z0 = zs(t)
    seeds = rng.integers(0, 2**31, NZ)
    zb = [x for x in Parallel(n_jobs=args.jobs)(delayed(C.safe)(lambda s: zs(C.boot(t, np.random.default_rng(s))), s) for s in seeds) if x is not None]
    zb = [z for z in zb if all(k in z for k in KC)]
    pd.DataFrame([{f"{k}_{p}": getattr(z[k], p) for k in KC for p in ("real", "imag")} for z in zb]).to_csv(OUT / "z_boot.csv", index=False)
    print("z0", {k: (abs(z0[k]).round(3), (np.angle(z0[k]) * HOURS).round(2)) for k in KC}, "n_boot", len(zb), flush=True)
    # ---- gains
    Lp = M.bulk_projector(M.bulk())
    pb_ps, _ = C.load("psoriasis_sc_bins")
    wps, _ = bin_weights(pb_ps)
    hb = pbh[(pbh.status == "healthy") & pbh.compartment.isin(KC) & (pbh.n_cells >= C.MIN_CELLS) & (pbh.stratum != "GSE162183")]
    xh = hb.groupby("compartment"); dep = xh.depth.sum().reindex(KC)
    cpm_h = xh[H.GENES].sum().reindex(KC).div(dep, axis=0).T * 1e6
    CPMK = {"lesional": wps["L"][1][KC], "nonlesional": wps["N"][1][KC], "healthy": cpm_h}
    def gains_of(means):   # {group: Series of shares} -> {group: gains}
        return {g: M.gains(M.shares(means[g], CPMK[g]), Lp) for g in GROUPS}
    steps = {"L-U": ("lesional", "nonlesional"), "U-H": ("nonlesional", "healthy")}
    def stats(means_by_series, keys, z):
        """-> {(step, series|pooled): offset h}"""
        out = {}
        for step, (a, b) in steps.items():
            rat = {}
            for s in keys:
                g = gains_of({gg: means_by_series[s][gg] for gg in GROUPS})
                rat[s] = ratio(g[a], g[b], z); out[(step, s)] = float(np.angle(rat[s]) * HOURS)
            out[(step, "pooled")] = pool_unit([rat[s] for s in keys], [w54.get(s, 1.0) for s in keys])
        return out
    results = []
    obs_tab = {}
    for tagset, keys in (("primary_3series", [s for s in SERIES if s != "GSE13355"]), ("sensitivity_4series", SERIES)):
        means = {s: {g: SH[(s, g)].mean() for g in GROUPS} for s in keys}
        obs = stats(means, keys, z0)
        rr = np.random.default_rng(SEED + (1 if tagset.startswith("sens") else 0))
        boot = []
        for i in range(NB):
            mb = {s: {g: SH[(s, g)].iloc[rr.integers(0, len(SH[(s, g)]), len(SH[(s, g)]))].mean() for g in GROUPS} for s in keys}
            boot.append(stats(mb, keys, zb[i % len(zb)]))
        B = pd.DataFrame(boot)
        for k in obs:
            bv = B[k].to_numpy(); v = obs[k]
            d = np.asarray(A.S.wrap(bv - v)); lo, hi = v + np.quantile(d, .025), v + np.quantile(d, .975)
            results.append(dict(set=tagset, step=k[0], series=k[1], predicted_h=v, ci_lo=lo, ci_hi=hi, boot_sd=bv.std(ddof=1), n_boot=NB,
                                p_boot_neg=(1 + np.sum(bv >= 0)) / (1 + NB), p_boot_pos=(1 + np.sum(bv <= 0)) / (1 + NB)))
        print(tagset, "done", flush=True)
        obs_tab[tagset] = observed(items, sets, w54, keys, tagset)
    res = pd.DataFrame(results)
    # ---- hypotheses (pooled)
    h = []
    for tagset in res.set.unique():
        r = res[(res.set == tagset) & (res.series == "pooled")].set_index("step")
        p = np.array([r.loc["L-U", "p_boot_neg"], r.loc["U-H", "p_boot_pos"]]); padj = holm(p)
        for i, (hid, step, sign) in enumerate((("H148.1", "L-U", "<0"), ("H148.2", "U-H", ">0"))):
            v, lo, hi = r.loc[step, ["predicted_h", "ci_lo", "ci_hi"]]
            hit = (v < 0 and hi < 0) if sign == "<0" else (v > 0 and lo > 0)
            h.append(dict(set=tagset, hypothesis=hid, step=step, direction=sign, predicted_h=v, ci_lo=lo, ci_hi=hi, p_one_sided=p[i],
                          p_holm=padj[i], bound_ok=hit, verdict="PASS" if hit and padj[i] < 0.05 else "FAIL",
                          acceptance_pass=ok[tagset]))
    hd = pd.DataFrame(h); hd.to_csv(OUT / "h148.csv", index=False)
    res.to_csv(OUT / "h148_predicted_all.csv", index=False)
    print(hd.round(4).to_string(), flush=True)
    # ---- D148
    d = []
    for tagset in res.set.unique():
        for step in ("L-U", "U-H"):
            for s in list(res[res.set == tagset].series.unique()):
                r = res[(res.set == tagset) & (res.step == step) & (res.series == s)].iloc[0]
                o = obs_tab[tagset][(step, s)]; pr = r.predicted_h
                d.append(dict(set=tagset, step=step, series=s, observed_h=o, predicted_h=pr, pred_ci_lo=r.ci_lo, pred_ci_hi=r.ci_hi,
                              share_explained=pr / o if abs(o) > 1e-9 else np.nan, residual_h=float(A.wrap(o - pr))))
    pd.DataFrame(d).to_csv(OUT / "d148.csv", index=False)
    print(pd.DataFrame(d).round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
