"""A62 and A77: layer-resolved psoriatic delay (epidermis vs dermis) and the basal layer in intact tissue.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A62 H62.1, A77 H77.1).
A62: layer-resolved spatial data (Visium GSE206391, GSE202011; GeoMx GSE310145), one pseudobulk per patient x status x
layer. The A54 estimator adapted to small series (fixed before scoring): one projector per series and layer, states
centred within series and status, the pooled complex correlation from cross-moments summed over series (weights =
pairs), the reference gate on the series-centred uninvolved epidermis of all series, one-sided per-pair swap null
(10,000). H62.1: epidermal offset < 0, p < 0.05 (UNINFORMATIVE if the gate fails). Secondary: dermis, epidermis minus
dermis, variance ratio; sensitivities: psoriasis vulgaris only, GeoMx epidermis segment only, without GeoMx, leave one
series out.
A77: GeoMx GSE310145, basal layer vs suprabasal epidermis vs dermis of the same patient x status (6 patients, lesional
sites summed). Offset of one compartment against another = arg of the complex correlation (h, > 0: the first ahead),
patient bootstrap; compression (lesional minus uninvolved offset) against all 2^6 patient-level status swaps. H77.1:
uninvolved epidermis minus basal > 0 with the bootstrap interval > 0.

    python compute_layers.py [--subsample N]
Outputs: paper/results/figure_5/output/psoriasis_details/a62_layers.csv, a62_series.csv, a62_verdict.csv, a77_offsets.csv,
  a77_compression.csv, a77_verdict.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a62_layers.py, a77_geomx_layers.py
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402

PD, S = H.PD, H.S
SERIES = ["GSE206391", "GSE202011", "GSE310145"]


def load(layers, gse, layer, **kw):
    """(case, ref) aligned on patients, or None (legacy `a62_layers.load`)."""
    t = layers.layer_pairs(SAMPLES, layer, series=[gse], **kw)
    if t.empty:
        return None
    case, ref = H.split_pairs(t)
    return (case, ref) if len(case) else None


def pooled_corr(cs) -> complex:
    """local: not in clocksync. Complex correlation from cross-moments summed over series (states centred within
    series), Poisson-corrected denominators."""
    c = [(S.complex_state(bl), S.complex_state(bn), vl, vn) for bl, vl, bn, vn in cs]
    num = sum(np.sum(a * np.conj(b)) for a, b, _, _ in c)
    el = sum(np.sum(abs(a) ** 2) - (len(a) - 1) * np.mean(vl) for a, _, vl, _ in c)
    en = sum(np.sum(abs(b) ** 2) - (len(b) - 1) * np.mean(vn) for _, b, _, vn in c)
    return num / np.sqrt(max(el, 1e-12) * max(en, 1e-12))


def swap(cs, rng):
    out = []
    for bl, vl, bn, vn in cs:
        s = rng.random(len(bl)) < 0.5
        out.append((np.where(s[:, None], bn, bl), np.where(s, vn, vl), np.where(s[:, None], bl, bn),
                    np.where(s, vl, vn)))
    return out


def analyse(layers, layer, label, n_perm, series=SERIES, **kw):
    items = []
    for g in series:
        x = load(layers, g, layer, **kw)
        if x is not None and len(x[0]) >= 2:
            items.append((g, *x))
    cs = [PD.cohort_states(c, r, True, H.WF) for _, c, r in items]
    c = pooled_corr(cs)
    off = float(np.angle(c) * H.HOURS)
    rng = np.random.default_rng(62)
    null = np.array([np.angle(pooled_corr(swap(cs, rng))) * H.HOURS for _ in range(n_perm)])
    w = np.array([len(c_) for _, c_, _ in items], float)
    r = float(np.sum(w * [PD.log_var_ratio(*b) for b in cs]) / w.sum())
    v = dict(layer=layer, set=label, series=";".join(g for g, _, _ in items), n_pairs=int(w.sum()), offset_h=off,
             corr_abs=float(abs(c)), p_lower=D.p_lower(null, off), r_mean=r)
    per = [dict(layer=layer, set=label, series=g, n_pairs=len(c_), offset_h=H.offset_h(b), r=PD.log_var_ratio(*b))
           for (g, c_, _), b in zip(items, cs)]
    return v, per


def gate(layers):
    ys = []
    for g in SERIES:
        x = load(layers, g, "epidermis")
        if x is None:
            continue
        y = H.logexpr(x[1], H.GENES, True)
        ys.append(y - y.mean(0))
    return PD.reference_gate(np.vstack(ys), H.WF, 200, np.random.default_rng(62))


def a62(layers, out, n_perm):
    gt = gate(layers)
    print("gate", gt, flush=True)
    specs = [("epidermis", "primary", {}), ("dermis", "primary", {}),
             ("epidermis", "psoriasis_vulgaris", dict(pso_only=True)),
             ("dermis", "psoriasis_vulgaris", dict(pso_only=True)),
             ("epidermis", "geomx_epidermis_segment", dict(regions=("epidermis",))),
             ("epidermis", "no_GeoMx", dict(series=SERIES[:2])), ("dermis", "no_GeoMx", dict(series=SERIES[:2]))]
    specs += [("epidermis", f"without_{g}", dict(series=[s for s in SERIES if s != g])) for g in SERIES]
    rows, per = [], []
    for layer, label, kw in specs:
        v, p = analyse(layers, layer, label, n_perm, **kw)
        rows.append(v)
        per += p
        print({k: (round(x, 3) if isinstance(x, float) else x) for k, x in v.items()}, flush=True)
    d = pd.DataFrame(rows)
    o = d.set_index(["layer", "set"]).offset_h
    v = dict(gate_share=gt["share"], gate_p=gt["p"], gate_pass=gt["p"] < 0.05,
             epidermis_minus_dermis_h=float((o[("epidermis", "primary")] - o[("dermis", "primary")] + 12) % 24 - 12))
    pe = d[(d.layer == "epidermis") & (d.set == "primary")].iloc[0]
    v["H62_1"] = ("UNINFORMATIVE (gate)" if not v["gate_pass"] else
                  "PASS" if (pe.offset_h < 0 and pe.p_lower < 0.05) else "FAIL")
    d.to_csv(out / "a62_layers.csv", index=False)
    pd.DataFrame(per).to_csv(out / "a62_series.csv", index=False)
    pd.Series(v).to_csv(out / "a62_verdict.csv", header=False)
    print(d.round(3).to_string(), "\n", pd.Series(v).to_string(), flush=True)


REGIONS = {"basal": ("basal_layer",), "epidermis": ("epidermis",), "dermis": ("vim_high", "vim_low")}
PAIRS = [("epidermis", "basal"), ("basal", "dermis"), ("epidermis", "dermis")]


def geomx_tables(samples) -> dict:
    d = samples[samples.series == "GSE310145"]
    out = {}
    for comp, regs in REGIONS.items():
        x = d[d.region.isin(regs)]
        g = x.groupby(["condition", "unit"])
        t = g[H.GENES].sum()
        t["depth"] = g.depth.sum()
        out[comp] = t
    return out


def offset(a, b, idx) -> float:
    """Offset of compartment a relative to b over units idx (one projector)."""
    return float(np.angle(PD.ccorr(*PD.cohort_states(a.loc[idx], b.loc[idx], True, H.WF))) * H.HOURS)


def a77(samples, out, n_boot):
    T = geomx_tables(samples)
    pts = sorted(set.intersection(*[set(t.xs(s, level="condition").index) for t in T.values()
                                    for s in ("lesional", "nonlesional")]))
    rng = np.random.default_rng(77)
    rows, comp = [], {}
    for a, b in PAIRS:
        for st in ("nonlesional", "lesional"):
            A, B = T[a].xs(st, level="condition"), T[b].xs(st, level="condition")
            o = offset(A, B, pts)
            bs = np.array([offset(A, B, list(rng.choice(pts, len(pts)))) for _ in range(n_boot)])
            dd = np.asarray((bs - o + 12) % 24 - 12)
            rows.append(dict(pair=f"{a}-{b}", status=st, n=len(pts), offset_h=o, lo=o + np.quantile(dd, 0.025),
                             hi=o + np.quantile(dd, 0.975)))

        def diff(flip):
            res = []
            for st in ("lesional", "nonlesional"):
                pa = {p: ("nonlesional" if st == "lesional" else "lesional") if p in flip else st for p in pts}
                A = pd.DataFrame([T[a].loc[(pa[p], p)] for p in pts], index=pts)
                B = pd.DataFrame([T[b].loc[(pa[p], p)] for p in pts], index=pts)
                res.append(offset(A, B, pts))
            return (res[0] - res[1] + 12) % 24 - 12

        obs = diff(set())
        null = np.array([diff(set(f)) for r in range(len(pts) + 1) for f in itertools.combinations(pts, r)])
        comp[f"{a}-{b}"] = dict(lesional_minus_uninvolved_h=obs, p_lower=float(np.mean(null <= obs)),
                                n_swaps=len(null))
    d, c = pd.DataFrame(rows), pd.DataFrame(comp).T
    e = d[(d.pair == "epidermis-basal") & (d.status == "nonlesional")].iloc[0]
    v = dict(H77_1="PASS" if (e.offset_h > 0 and e.lo > 0) else "FAIL")
    d.to_csv(out / "a77_offsets.csv", index=False)
    c.to_csv(out / "a77_compression.csv")
    pd.Series(v).to_csv(out / "a77_verdict.csv", header=False)
    print(d.round(2).to_string(), "\n", c.round(3).to_string(), "\n", pd.Series(v).to_string(), flush=True)


def main(argv=None):
    global SAMPLES
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    layers = H.builder("skin_layers")
    samples, _ = layers.build()
    SAMPLES = H.thin(samples, args.subsample, args.seed)
    sz = H.sizes(args, run, 10000, 1000)
    a62(layers, out, sz["n_perm"])
    a77(samples, out, sz["n_boot"])


if __name__ == "__main__":
    main()
