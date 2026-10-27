"""A115 scoring (registered 2026-10-27): trait component of the psoriatic delay, descriptive.

Per series one projector for the uninvolved and healthy rows together (`_lib.states_for`, GTEx core waveform and GTEx cue set A81b);
group offset = arg E[c_core c_cue*] over the group's samples (`S.ccorr`); pooled offset per group = arg of the A54-weighted sum of
the series' complex group statistics; D115.1 = uninvolved minus healthy (h, wrapped). 95 % bootstrap (2,000; samples resampled
within group and series, projector fixed), one-sided permutation p (10,000; group labels permuted within series; p = (1 + #{null
<= observed}) / (1 + n)), coherence per group (|weighted mean complex correlation|), atopic dermatitis (GSE121212), keratinocyte
differentiation score (A61 sets: DIFF minus BASAL, z-scored genes then score within series over uninvolved + healthy) as the
uninvolved-minus-healthy standardised mean difference (Hedges g; pooled = inverse-variance fixed effect). Seed 115.

    python score_a115.py     (needs output/a115_acceptance.csv and the parquet of build_healthy.py)
Outputs: output/a115_d115_1.csv, a115_atopic.csv, a115_keratinocyte.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import _lib as L
from clocksync.config import limit_threads

limit_threads(4)
OUT = Path(__file__).resolve().parents[1] / "output"
ND, NB, NP, SEED = None, 2000, 10000, 115
HOURS = L.HOURS
BASAL = "KRT5 KRT14 KRT15 COL17A1 TP63 ITGA6 ITGB4 CXCL14".split()
DIFF = "KRT1 KRT10 KRT2 DSG1 DSC1 SBSN KRTDAP LORICRIN FLG IVL CNFN LCE1B LCE3D SPRR1A SPRR1B SPRR2A".split()
ALIAS = {"LORICRIN": "LOR"}
sets = L.cue_sets(); w54 = L.weights54()
genes = sorted(set(L.all_set_genes(sets)) | set(DIFF) | set(BASAL) | set(ALIAS.values()))
DIR = L.ROOT / "processed_a115"


def prep(s, dis):
    d = L.read_series(DIR / f"{s}_{dis}.parquet", genes)
    counts = bool((d.value_type == "counts").all())
    x = pd.concat([d[d.status == "nonlesional"].drop_duplicates("patient"), d[d.status == "healthy"].drop_duplicates("patient")],
                  ignore_index=True)
    grp = (x.status == "healthy").to_numpy()          # True = healthy
    bc, vc = L.states_for(x, L.WF, counts)[0]; bk, vk = L.states_for(x, sets["cue"], counts)[0]
    from clocksync import paired as PD
    names = {g: g for g in x.columns}
    def sc(pos, neg):
        lx = PD.log_expr(x[[g for g in x.columns if g not in L.META]].astype(float), x.depth.to_numpy(float) if counts else None)
        z = (lx - lx.mean()) / lx.std().replace(0, 1)
        f = lambda gs: [g if g in z else ALIAS.get(g, g) for g in gs]
        s_ = z[[g for g in f(pos) if g in z]].mean(1) - z[[g for g in f(neg) if g in z]].mean(1)
        return ((s_ - s_.mean()) / s_.std()).to_numpy()
    return dict(n_u=int((~grp).sum()), n_h=int(grp.sum()), grp=grp, st=(bc, vc, bk, vk), kc=sc(DIFF, BASAL),
                kc_n=(len([g for g in DIFF if g in x or ALIAS.get(g) in x]), len([g for g in BASAL if g in x])))


def cc(st, idx):
    bc, vc, bk, vk = st
    return L.S.ccorr(bc[idx], bk[idx], vc[idx], vk[idx])


def stat(items, idxs, keys):
    """idxs[key] = (idx_U, idx_H) -> (pooled U offset, H offset, D, cU, cH) with A54 weights."""
    w = np.array([w54.get(k.split("|")[0], 1.0) for k in keys])
    cU = np.array([cc(items[k]["st"], idxs[k][0]) for k in keys]); cH = np.array([cc(items[k]["st"], idxs[k][1]) for k in keys])
    U, H = np.sum(w * cU), np.sum(w * cH)
    return L.wrap((np.angle(U) - np.angle(H)) * HOURS), np.angle(U) * HOURS, np.angle(H) * HOURS, abs(U) / w.sum(), abs(H) / w.sum()


def analyse(keys, items, tag, rng):
    base = {k: (np.where(~items[k]["grp"])[0], np.where(items[k]["grp"])[0]) for k in keys}
    obs = stat(items, base, keys)
    boot = np.array([stat(items, {k: (rng.choice(base[k][0], len(base[k][0])), rng.choice(base[k][1], len(base[k][1]))) for k in keys}, keys)[0]
                     for _ in range(NB)])
    null = np.empty(NP)
    for i in range(NP):
        idx = {}
        for k in keys:
            p = rng.permutation(len(items[k]["grp"])); nu = len(base[k][0]); idx[k] = (p[:nu], p[nu:])
        null[i] = stat(items, idx, keys)[0]
    lo, hi = np.quantile(boot, [0.025, 0.975])
    row = dict(analysis=tag, series="pooled", n_uninvolved=sum(items[k]["n_u"] for k in keys), n_healthy=sum(items[k]["n_h"] for k in keys),
               D115_h=obs[0], ci_lo=lo, ci_hi=hi, boot_n=NB, boot_sd=boot.std(ddof=1), p_one_sided_le=(1 + np.sum(null <= obs[0])) / (1 + NP),
               n_perm=NP, perm_null_sd=null.std(ddof=1), uninvolved_offset_h=obs[1], healthy_offset_h=obs[2], coherence_uninvolved=obs[3],
               coherence_healthy=obs[4], mde_h=(1.645 + 0.842) * null.std(ddof=1), n_series=len(keys))
    rows = [row]
    for k in keys:   # per series: same statistics, one series
        o1 = stat(items, base, [k]); sub = {k: items[k]}
        b1 = np.array([stat(items, {k: (rng.choice(base[k][0], len(base[k][0])), rng.choice(base[k][1], len(base[k][1])))}, [k])[0] for _ in range(NB)])
        n1 = np.empty(NP)
        for i in range(NP):
            p = rng.permutation(len(items[k]["grp"])); nu = len(base[k][0]); n1[i] = stat(items, {k: (p[:nu], p[nu:])}, [k])[0]
        l1, h1 = np.quantile(b1, [0.025, 0.975])
        rows.append(dict(analysis=tag, series=k, n_uninvolved=items[k]["n_u"], n_healthy=items[k]["n_h"], D115_h=o1[0], ci_lo=l1, ci_hi=h1,
                         boot_n=NB, boot_sd=b1.std(ddof=1), p_one_sided_le=(1 + np.sum(n1 <= o1[0])) / (1 + NP), n_perm=NP,
                         perm_null_sd=n1.std(ddof=1), uninvolved_offset_h=o1[1], healthy_offset_h=o1[2], coherence_uninvolved=o1[3],
                         coherence_healthy=o1[4], mde_h=(1.645 + 0.842) * n1.std(ddof=1), n_series=1,
                         weight_A54=w54.get(k.split("|")[0], np.nan)))
    return pd.DataFrame(rows)


def kerat(keys, items, tag):
    rows = []
    for k in keys:
        x, g = items[k]["kc"], items[k]["grp"]
        u, h = x[~g], x[g]; nu, nh = len(u), len(h)
        sp = np.sqrt(((nu - 1) * u.var(ddof=1) + (nh - 1) * h.var(ddof=1)) / (nu + nh - 2))
        J = 1 - 3 / (4 * (nu + nh) - 9); d = (u.mean() - h.mean()) / sp; gg = J * d
        se = J * np.sqrt((nu + nh) / (nu * nh) + d ** 2 / (2 * (nu + nh)))
        rows.append(dict(analysis=tag, series=k, n_uninvolved=nu, n_healthy=nh, smd_hedges_g=gg, se=se, ci_lo=gg - 1.96 * se, ci_hi=gg + 1.96 * se,
                         n_genes_diff=items[k]["kc_n"][0], n_genes_basal=items[k]["kc_n"][1]))
    d = pd.DataFrame(rows)
    if len(d) > 1:
        w = 1 / d.se ** 2; m = np.sum(w * d.smd_hedges_g) / w.sum(); se = 1 / np.sqrt(w.sum())
        d = pd.concat([pd.DataFrame([dict(analysis=tag, series="pooled", n_uninvolved=d.n_uninvolved.sum(), n_healthy=d.n_healthy.sum(),
                                          smd_hedges_g=m, se=se, ci_lo=m - 1.96 * se, ci_hi=m + 1.96 * se)]), d], ignore_index=True)
    return d


def main():
    acc = pd.read_csv(OUT / "a115_acceptance.csv")
    ok = acc.groupby(["series", "disease"]).accepted.all()
    dropped = [f"{s}|{d}" for (s, d), v in ok.items() if not v]
    print("dropped (failed acceptance):", dropped, flush=True)
    rng = np.random.default_rng(SEED)
    pso = [s for s in ("GSE13355", "GSE121212", "GSE201827", "GSE268958") if f"{s}|psoriasis" not in dropped]
    items = {s: prep(s, "psoriasis") for s in pso}
    for s in pso: print(s, items[s]["n_u"], items[s]["n_h"], flush=True)
    r = analyse(pso, items, "psoriasis", rng); r.to_csv(OUT / "a115_d115_1.csv", index=False); print(r.round(3).to_string(), flush=True)
    ad = {"GSE121212": prep("GSE121212", "atopic_dermatitis")}
    ra = analyse(["GSE121212"], ad, "atopic_dermatitis", rng); ra.to_csv(OUT / "a115_atopic.csv", index=False); print(ra.round(3).to_string(), flush=True)
    kc = pd.concat([kerat(pso, items, "psoriasis"), kerat(["GSE121212"], ad, "atopic_dermatitis")], ignore_index=True)
    kc.to_csv(OUT / "a115_keratinocyte.csv", index=False); print(kc.round(3).to_string())


if __name__ == "__main__":
    main()
