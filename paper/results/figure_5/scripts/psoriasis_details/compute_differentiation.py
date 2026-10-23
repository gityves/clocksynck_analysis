"""A63-A66: the keratinocyte clock along differentiation (healthy single-cell skin), its compression in psoriatic lesions,
and its dependence on body site.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A63 H63.1; A64, A66 descriptive; A65 H65.1, H65.2). Pseudobulks per donor x
compartment: endothelial, fibroblast, mural, kc_cycling and non-cycling keratinocytes in within-dataset quartiles of
(differentiated - basal score), kc_d1 (basal) ... kc_d4 (builders `healthy_skin_bins`, `psoriasis_sc_bins`,
`healthy_skin_sites`); >= 20 cells; free host (`Dataset(anchors=None)`), offsets relative to kc_d1 or to the dermal classes.
  A63 (healthy donors; GSE162183 excluded, no PER1): H63.1 offset(kc_d4) - offset(kc_d1) > 0, one-sided p against bin labels
    permuted among kc_d1..d4 within each donor (n_perm 1,000); secondary Spearman of bin rank vs offset, kc_cycling - kc_d1,
    donor bootstrap (500); uninvolved psoriatic skin ('nonlesional') and all donors descriptive. Gate per set.
  A64 (descriptive): offsets of kc_cycling, kc_d1..d4 relative to the circular mean of fibroblast, endothelial, mural and to
    fibroblast; share of the kc_d4-to-dermis offset already present at kc_d1; layer-balanced fit (dermal classes vs kc_d1,
    kc_d4, kc_cycling), donor bootstrap (500).
  A65 (psoriasis units with bins in both statuses): G = offset(kc_d4) - offset(kc_d1), E = circular mean of kc_d1..d4 minus
    circular mean of the dermal classes; H65.1 G_L - G_N < 0, H65.2 E_L - E_N < 0, one-sided p against statuses swapped per
    unit (1,000); the gate must pass; sensitivity without GSE278330; unit bootstrap (500).
  A66 (descriptive, healthy skin by body site; unit = donor x site): separate fits for volar (palm + sole), non-volar (hip +
    face), face, hip; offsets relative to fibroblast, G and E, classes in >= 3 units, unit bootstrap (300).

    python compute_differentiation.py [--subsample N] [--n-perm P] [--n-boot B] [--jobs 4] [--only a63|a64|a65|a66]
Outputs (paper/results/figure_5/output/psoriasis_details/): a63_{healthy,nonlesional,all}_compartments.csv, a63_verdict.csv,
a64_{offsets,balanced_fit}.csv, a65_{all,no_GSE278330}_offsets.csv, a65_verdict.csv, a66_sites.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a63_differentiation.py, a64_epi_dermis.py, a65_compression.py,
a66_sites.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import _details as D  # noqa: E402
import _shared as H  # noqa: E402
import _singlecell as C  # noqa: E402
from common import cli  # noqa: E402

DERMAL, BINS = C.DERMAL, C.BINS
COMP = ["endothelial", "fibroblast", "mural", "kc_cycling"] + BINS
KC = ["kc_cycling"] + BINS
BALANCED = DERMAL + ["kc_cycling", "kc_d1", "kc_d4"]


# ---- tables ---------------------------------------------------------------------------------------------------------

def healthy_table(pb, status):
    """A63 table: donors of the given statuses, one unit per study | patient | visit | status, >= 5 units per compartment."""
    d = pb[pb.status.isin(status) & pb.compartment.isin(COMP) & (pb.n_cells >= C.MIN_CELLS) & (pb.stratum != "GSE162183")]
    n = d.groupby("compartment").donor.nunique()
    return d[d.compartment.isin(n.index[n >= C.MIN_UNITS])].reset_index(drop=True)


def permuted(t, rng):
    """Bin labels permuted among kc_d1..d4 within each donor (A63 null)."""
    s = t.copy()
    for _, idx in s[s.compartment.isin(BINS)].groupby("donor").groups.items():
        s.loc[idx, "compartment"] = rng.permutation(s.loc[idx, "compartment"].to_numpy())
    return s


# ---- A63 ------------------------------------------------------------------------------------------------------------

def stat63(t) -> dict:
    tab = C.fit_states(t)
    off = {k: float(C.wrap(tab.offset_h[k] - tab.offset_h["kc_d1"])) for k in BINS + ["kc_cycling"] if k in tab.index}
    bins = [k for k in BINS if k in off]
    rho = spearmanr(range(len(bins)), [off[k] for k in bins]).statistic if len(bins) >= 3 else np.nan
    return dict(d4_minus_d1=off.get("kc_d4", np.nan), cycling_minus_d1=off.get("kc_cycling", np.nan), rho=rho, rel=off,
                tab=tab)


def _s63(f, t, seed):
    out = stat63(f(t, np.random.default_rng(seed)))
    out.pop("tab")
    return out


def run63(t, tag, n, jobs, out, test=False, seed=63):
    gate = C.gate(t, seed)
    obs = stat63(t)
    seeds = np.random.default_rng(seed).integers(0, 2**31, n["n_perm"] + n["n_boot"])
    par = Parallel(n_jobs=jobs)
    boot = [b for b in par(delayed(C.safe)(_s63, C.boot, t, s) for s in seeds[n["n_perm"]:]) if b]
    tab = obs.pop("tab")
    rows = []
    for k in tab.index:
        r = dict(compartment=k, offset_h=tab.offset_h[k], coupling=tab.coupling[k], amplitude=tab.rel_amplitude[k],
                 n_units=int(t[t.compartment == k].donor.nunique()))
        if k in obs["rel"]:
            lo, hi = C.interval([b["rel"].get(k, np.nan) for b in boot], obs["rel"][k])
            r.update(rel_d1_h=obs["rel"][k], rel_lo=lo, rel_hi=hi)
        rows.append(r)
    v = dict(set=tag, donors=t.donor.nunique(), datasets=t.stratum.nunique(), gate_share=gate["clock_share"],
             gate_p_acrophase=gate["p_acrophase"], gate_p_pairing=gate["p_pairing"], gate_pass=gate["gate_pass"],
             d4_minus_d1_h=obs["d4_minus_d1"], cycling_minus_d1_h=obs["cycling_minus_d1"], spearman_rank=obs["rho"])
    for key in ("d4_minus_d1", "cycling_minus_d1"):
        v[f"{key}_lo"], v[f"{key}_hi"] = C.interval([b[key] for b in boot], obs[key])
    if test:
        null = [x for x in par(delayed(C.safe)(_s63, permuted, t, s) for s in seeds[:n["n_perm"]]) if x]
        null = np.array([x["d4_minus_d1"] for x in null], float)
        v["p_upper"] = float((1 + np.sum(null >= obs["d4_minus_d1"])) / (1 + np.isfinite(null).sum()))
        v["H63_1"] = ("UNINFORMATIVE (gate)" if not v["gate_pass"] else
                      "PASS" if (obs["d4_minus_d1"] > 0 and v["p_upper"] < 0.05) else "FAIL")
    pd.DataFrame(rows).to_csv(out / f"a63_{tag}_compartments.csv", index=False)
    print(tag, "\n", pd.DataFrame(rows).round(2).to_string(), "\n", pd.Series(v).to_string(), flush=True)
    return v


# ---- A64 ------------------------------------------------------------------------------------------------------------

def rel64(t) -> dict:
    off = C.fit_states(t).offset_h
    derm = C.cmean([off[k] for k in DERMAL if k in off.index])
    out = {f"{k}_vs_dermal": float(C.wrap(off[k] - derm)) for k in KC if k in off.index}
    out |= {f"{k}_vs_fibroblast": float(C.wrap(off[k] - off["fibroblast"])) for k in KC if k in off.index}
    return out


def run64(t, tag, n, jobs, seed=64):
    obs = rel64(t)
    seeds = np.random.default_rng(seed).integers(0, 2**31, n["n_boot"])
    bs = [b for b in Parallel(n_jobs=jobs)(delayed(C.safe)(lambda s: rel64(C.boot(t, np.random.default_rng(s))), s)
                                           for s in seeds) if b]
    rows = []
    for k, v in obs.items():
        lo, hi = C.interval([b.get(k, np.nan) for b in bs], v)
        rows.append(dict(set=tag, quantity=k, h=v, lo=lo, hi=hi))
    sh = [b["kc_d1_vs_dermal"] / b["kc_d4_vs_dermal"] for b in bs if abs(b.get("kc_d4_vs_dermal", 0)) > 0.1]
    rows.append(dict(set=tag, quantity="share_d1_of_d4", h=obs["kc_d1_vs_dermal"] / obs["kc_d4_vs_dermal"],
                     lo=np.quantile(sh, 0.025), hi=np.quantile(sh, 0.975)))
    bal = C.fit_states(t[t.compartment.isin(BALANCED)].reset_index(drop=True))
    return pd.DataFrame(rows), bal.assign(set=tag).reset_index()


# ---- A65 ------------------------------------------------------------------------------------------------------------

def stat65(t) -> dict:
    o = C.fit_states(t).offset_h
    out = {}
    for s in "LN":
        bins = [o[f"{s}:{k}"] for k in BINS if f"{s}:{k}" in o]
        out[f"G_{s}"] = float(C.wrap(o[f"{s}:kc_d4"] - o[f"{s}:kc_d1"]))
        out[f"E_{s}"] = float(C.wrap(C.cmean(bins) - C.cmean([o[f"{s}:{k}"] for k in DERMAL if f"{s}:{k}" in o])))
    out["dG"] = float(C.wrap(out["G_L"] - out["G_N"]))
    out["dE"] = float(C.wrap(out["E_L"] - out["E_N"]))
    return out


def run65(t, tag, n, jobs, out, test, seed=65):
    gate = C.gate(t, seed)
    obs = stat65(t)
    seeds = np.random.default_rng(seed).integers(0, 2**31, n["n_perm"] + n["n_boot"])
    par = Parallel(n_jobs=jobs)
    bs = [b for b in par(delayed(C.safe)(lambda s: stat65(C.boot(t, np.random.default_rng(s))), s)
                         for s in seeds[n["n_perm"]:]) if b]
    keys = ("G_L", "G_N", "dG", "E_L", "E_N", "dE")
    v = dict(set=tag, units=t.donor.nunique(), gate_share=gate["clock_share"], gate_p_acrophase=gate["p_acrophase"],
             gate_p_pairing=gate["p_pairing"], gate_pass=gate["gate_pass"], **{k: obs[k] for k in keys})
    for k in keys:
        v[f"{k}_lo"], v[f"{k}_hi"] = C.interval([b[k] for b in bs], obs[k])
    if test:
        null = [x for x in par(delayed(C.safe)(lambda s: stat65(C.swapped(t, np.random.default_rng(s))), s)
                               for s in seeds[:n["n_perm"]]) if x]
        for k in ("dG", "dE"):
            nv = np.array([x[k] for x in null], float)
            v[f"p_{k}"] = float((1 + np.sum(nv <= obs[k])) / (1 + len(nv)))
        ok = v["gate_pass"]
        v["H65_1"] = "UNINFORMATIVE (gate)" if not ok else "PASS" if obs["dG"] < 0 and v["p_dG"] < 0.05 else "FAIL"
        v["H65_2"] = "UNINFORMATIVE (gate)" if not ok else "PASS" if obs["dE"] < 0 and v["p_dE"] < 0.05 else "FAIL"
    offs = C.fit_states(t).offset_h.rename("offset_h").rename_axis("compartment").reset_index().assign(set=tag)
    offs.to_csv(out / f"a65_{tag}_offsets.csv", index=False)
    print(tag, "\n", offs.round(2).to_string(), "\n", pd.Series(v).to_string(), flush=True)
    return v


def table65(pb):
    d = C.status_split(pb, COMP)
    n = d.groupby("compartment").donor.nunique()
    ok = [c for c in COMP if n.get(f"L:{c}", 0) >= C.MIN_UNITS and n.get(f"N:{c}", 0) >= C.MIN_UNITS]
    return d[d.compartment.str[2:].isin(ok)].reset_index(drop=True)


# ---- A66 ------------------------------------------------------------------------------------------------------------

def table66(pb, donors, sel, min_units=3):
    d = pb.merge(donors[["site", "volar"]], left_on="donor", right_index=True)
    d["study"] = d.stratum
    d = d[sel(d) & d.compartment.isin(COMP) & (d.n_cells >= 20)]
    n = d.groupby("compartment").donor.nunique()
    return d[d.compartment.isin(n.index[n >= min_units])].reset_index(drop=True)


def stat66(t) -> dict:
    o = C.fit_states(t, min_donors=3).offset_h
    bins = [o[k] for k in BINS if k in o]
    out = {f"{k}_vs_fibroblast": float(C.wrap(o[k] - o["fibroblast"])) for k in KC if k in o}
    out["G"] = float(C.wrap(o["kc_d4"] - o["kc_d1"])) if {"kc_d1", "kc_d4"} <= set(o.index) else np.nan
    out["E"] = float(C.wrap(C.cmean(bins) - C.cmean([o[k] for k in DERMAL if k in o])))
    return out


GROUPS66 = {"volar": lambda d: d.volar.astype(bool), "non_volar": lambda d: ~d.volar.astype(bool),
            "face": lambda d: d.study.astype(str).str.contains("13085"),
            "hip": lambda d: d.site.astype(str).str.lower().eq("hip")}


def run66(args, n, out, seed=66):
    pb, donors = C.load("healthy_skin_sites", args.subsample, args.seed)
    rows = []
    for g, sel in GROUPS66.items():
        t = table66(pb, donors, sel)
        try:
            obs = stat66(t)
        except Exception as e:  # too sparse to fit
            rows.append(dict(group=g, units=t.donor.nunique(), note=f"not estimable: {e}"))
            continue
        seeds = np.random.default_rng(seed).integers(0, 2**31, n["n_boot"])
        bs = [b for b in Parallel(n_jobs=args.jobs)(delayed(C.safe)(lambda s: stat66(C.boot(t, np.random.default_rng(s))),
                                                                    s) for s in seeds) if b]
        for k, v in obs.items():
            bv = np.array([b.get(k, np.nan) for b in bs], float)
            lo, hi = C.interval(bv, v)
            rows.append(dict(group=g, units=t.donor.nunique(), donors=t.donor.str.split("|").str[1].nunique(),
                             compartments=";".join(sorted(t.compartment.unique())), quantity=k, h=v, lo=lo, hi=hi,
                             n_boot=int(np.isfinite(bv).sum())))
        print(pd.DataFrame([r for r in rows if r["group"] == g]).round(2).to_string(), flush=True)
    pd.DataFrame(rows).to_csv(out / "a66_sites.csv", index=False)


def main(argv=None):
    ap = H.parser(__doc__)
    ap.add_argument("--only", choices=("a63", "a64", "a65", "a66"), default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    n = dict(n_perm=args.n_perm or (30 if run["quick"] else 1000), n_boot=args.n_boot or (20 if run["quick"] else 500))
    pick = lambda a: args.only in (None, a)
    if pick("a63") or pick("a64"):
        pb, _ = C.load("healthy_skin_bins", args.subsample, args.seed)
        if pick("a63"):
            v = [run63(healthy_table(pb, ["healthy"]), "healthy", n, args.jobs, out, test=True),
                 run63(healthy_table(pb, ["nonlesional"]), "nonlesional", n, args.jobs, out),
                 run63(healthy_table(pb, ["healthy", "nonlesional"]), "all", n, args.jobs, out)]
            pd.DataFrame(v).to_csv(out / "a63_verdict.csv", index=False)
        if pick("a64"):
            r, b = zip(*[run64(healthy_table(pb, st), tag, n, args.jobs)
                         for st, tag in ((["healthy"], "healthy"), (["healthy", "nonlesional"], "all"))])
            pd.concat(r).to_csv(out / "a64_offsets.csv", index=False)
            pd.concat(b).to_csv(out / "a64_balanced_fit.csv", index=False)
    if pick("a65"):
        pb, _ = C.load("psoriasis_sc_bins", args.subsample, args.seed)
        t = table65(pb)
        v = [run65(t, "all", n, args.jobs, out, True),
             run65(t[t.stratum != "GSE278330"].reset_index(drop=True), "no_GSE278330", n, args.jobs, out, False)]
        pd.DataFrame(v).to_csv(out / "a65_verdict.csv", index=False)
    if pick("a66"):
        run66(args, dict(n, n_boot=args.n_boot or (20 if run["quick"] else 300)), out)


if __name__ == "__main__":
    main()
