"""A185 (C1, secondary): is A175's lower tumour-normal pairing in proliferating tumours explained by less tumour clock signal?

Registered in paper/docs/registered_experiments.md, section "A185". Reuses A175 (`run_a175`: setup_data, prepare, tertiles,
Contrast, subset, pseudo-set construction) and `q3_lib` / `h138_8.Cancer` by import.
Amplitude a = |z_tumour| (centred tumour state `st[c]["za"]`, the state entering `q3_lib.corr_stat`). Within each A175 cancer
(>= 10 pairs per arm), the low- and high-G2M-tertile pairs are pooled and cut into amplitude tertiles (rank). In each
cancer x amplitude bin with >= 5 pairs per arm, pairing(low) - pairing(high) (H138.4 statistic, states re-centred within the
subset); D_match = mean over bins weighted by n_low + n_high (common weights, so the arms are compared on the same bins).
Null: arm labels permuted within cancer x amplitude bin (10,000). Power stage uses labels only.
Companions: D_match on the 200 matched pseudo-clock and 200 responsive sets of A175 (same streams, same tertile labels);
amplitude bins re-cut on each set's own tumour amplitude.

    python run_a185.py --stage reproduce|power|tests [--jobs 4] [--n-perm 10000] [--n-sets 200]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parent / "A175_c1_pairing_strata" / "scripts"))
import run_a175 as r75  # noqa: E402  (also puts A138 scripts on the path)
import q3_lib as q  # noqa: E402
from h138_8 import Cancer  # noqa: E402
from common import cli  # noqa: E402
from power import SEED  # noqa: E402

MIN_BIN, NB = 5, 3
log = r75.log


def amp_of(C: Cancer) -> np.ndarray:
    """|centred tumour state| (= |za| of q3_lib.tcga_states)."""
    return np.abs(C.a - C.a.mean())


def amp_bins(lab: np.ndarray, amp: np.ndarray) -> np.ndarray:
    """-1 outside the low (0) / high (2) arms, else amplitude tertile (rank, pooled low+high pairs of the cancer)."""
    b = -np.ones(len(lab), int)
    sel = (lab == 0) | (lab == 2)
    b[sel] = r75.tertiles(amp[sel])
    return b


def bin_stats(C: Cancer, lab: np.ndarray, amp: np.ndarray):
    """Per amplitude bin: (n_low, n_high, pairing low, pairing high) on the true labels."""
    bn = amp_bins(lab, amp)
    out = []
    for k in range(NB):
        lo, hi = (bn == k) & (lab == 0), (bn == k) & (lab == 2)
        out.append((int(lo.sum()), int(hi.sum()), float(C.stat(lo[None, :])[0]) if lo.sum() > 1 else np.nan,
                    float(C.stat(hi[None, :])[0]) if hi.sum() > 1 else np.nan))
    return out


def d_match(rows: list) -> tuple:
    """rows: (n_low, n_high, s_low, s_high) over cancer x bin; bins with >= MIN_BIN per arm; returns (D_match, n bins, n pairs)."""
    r = [x for x in rows if x[0] >= MIN_BIN and x[1] >= MIN_BIN]
    if not r:
        return np.nan, 0, 0
    d = np.array([x[2] - x[3] for x in r]); w = np.array([x[0] + x[1] for x in r], float)
    ok = np.isfinite(d)
    return (float(np.sum(d[ok] * w[ok]) / w[ok].sum()) if ok.any() else np.nan), len(r), int(w.sum())


class Match:
    """Spread-matched contrast on the A175 H175.1 cancers; `amp` = amplitude per cancer (full-length vector)."""

    def __init__(self, cs, d, cidx, use, amp, tag):
        self.cs, self.cidx, self.use, self.tag = cs, cidx, use, tag
        self.lab = {c: d[c]["tert"] for c in use}
        self.amp = {c: amp[c] for c in use}
        self.bin = {c: amp_bins(self.lab[c], self.amp[c]) for c in use}
        self.cells = []  # (cancer, bin, idx, n_low, n_high) for eligible bins
        self.dropped = []
        for c in use:
            for k in range(NB):
                idx = np.where(self.bin[c] == k)[0]
                nl, nh = int((self.lab[c][idx] == 0).sum()), int((self.lab[c][idx] == 2).sum())
                (self.cells if min(nl, nh) >= MIN_BIN else self.dropped).append((c, k, idx, nl, nh))

    def null(self, n_perm):
        rows = []
        for j, (c, k, idx, nl, nh) in enumerate(self.cells):
            rng = np.random.default_rng([SEED, 185, self.tag, self.cidx[c], k])
            labs = self.lab[c][idx]
            lp = labs[np.argsort(rng.random((n_perm, len(idx))), axis=1)]
            n = self.cs[c].n
            M = np.zeros((n_perm, n), bool); M[:, idx] = lp == 0
            H = np.zeros((n_perm, n), bool); H[:, idx] = lp == 2
            rows.append((self.cs[c].stat(M) - self.cs[c].stat(H), nl + nh))
        return q.pool(np.array([r[0] for r in rows]), [r[1] for r in rows])

    def observed(self):
        rows = []
        for c in self.use:
            rows += bin_stats(self.cs[c], self.lab[c], self.amp[c])
        return d_match(rows), rows


def amp_dict(cs, st, which):
    return {c: np.abs(st[c][which]) for c in cs}


def companion(c, i, n_sets, lab):
    """Per set and kind: bin rows (n_low, n_high, s_low, s_high) of cancer c, bins cut on the set's own tumour amplitude."""
    from clocksync.config import limit_threads
    from clocksync.paired import cohort_states, pseudo_waveforms
    limit_threads(1)
    samples = q.sh.load_tcga(None, 0)
    case, _ = q.sh.pair_frames(samples, c)
    ct, cn = q.sh.tcga_full_counts(c, case.index)
    ct, cn = ct.loc[case.index], cn.loc[case.index]
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    out = {}
    for kind, tag, kw in (("matched", 44, {}), ("responsive", 45, {"kind": "responsive"})):
        ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, tag, 0, i]), q.WF, **kw)
        res = []
        for w in ws:
            g = list(w.genes)
            bt, vt, bn, vn = cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)
            C = Cancer(c, dict(n=len(case), raw=(bt, vt, bn, vn)))
            res.append(bin_stats(C, lab, amp_of(C)))
        out[kind] = res
    return out


def setup():
    st, keep_c, cs, cidx, rep = r75.setup_data()
    d = r75.prepare(st, keep_c)
    con = r75.Contrast(cs, {c: d[c]["tert"] for c in cs}, 0, 2, 1, cidx)
    for c in keep_c:  # amplitude definition check
        assert np.allclose(amp_of(cs[c]), np.abs(st[c]["za"]), atol=1e-9), "amplitude != |za|"
    return st, keep_c, cs, cidx, rep, d, con


def stage_reproduce(args, out):
    st, keep_c, cs, cidx, rep, d, con = setup()
    _, pa, pb, D = con.observed()
    log(f"A175 H175.1 D reproduced: {D:.6f} (low {pa:.4f}, high {pb:.4f})")
    assert abs(D - 0.15773955999955447) < 1e-6, "A175 D not reproduced"
    pd.DataFrame([dict(D=D, pooled_low=pa, pooled_high=pb, a175_D=0.15773955999955447)]).to_csv(out / "reproduction_a175.csv", index=False)


def stage_power(args, out):
    st, keep_c, cs, cidx, rep, d, con = setup()
    stage_reproduce(args, out)
    m = Match(cs, d, cidx, con.use, amp_dict(cs, st, "za"), 1)
    null = m.null(args.n_perm)
    s = q.summarise(null)
    rows = [dict(hypothesis="H185.1", **{k: s[k] for k in ("null_mean", "null_q95", "mde", "n_perm")}, null_sd=s["null_sd"],
                 n_bins=len(m.cells), n_bins_dropped=len(m.dropped), n_pairs_in_bins=sum(x[3] + x[4] for x in m.cells),
                 not_estimable=bool(s["mde"] > 0.3), reference_D_A175=0.1577)]
    pd.DataFrame(rows).to_csv(out / "power.csv", index=False)
    pd.DataFrame([dict(cancer=c, bin=k, n_low=nl, n_high=nh, kept=True) for c, k, _, nl, nh in m.cells] +
                 [dict(cancer=c, bin=k, n_low=nl, n_high=nh, kept=False) for c, k, _, nl, nh in m.dropped]).to_csv(out / "bins.csv", index=False)
    log("power", rows[0])


def stage_tests(args, out):
    pw = pd.read_csv(out / "power.csv").iloc[0]
    st, keep_c, cs, cidx, rep, d, con = setup()
    _, pa, pb, D = con.observed()
    m = Match(cs, d, cidx, con.use, amp_dict(cs, st, "za"), 1)
    null = m.null(args.n_perm)
    s = q.summarise(null)
    assert np.isclose(s["mde"], pw.mde, rtol=1e-9), "null differs from power.csv"
    (Dm, nb, npairs), rows = m.observed()
    x = null[np.isfinite(null)]
    p = float((1 + np.sum(x >= Dm)) / (len(x) + 1))
    log(f"D_match {Dm:.4f} p {p:.4f} MDE {s['mde']:.4f} (D {D:.4f})")
    # companions
    cl = {c: i for i, c in enumerate(keep_c)}
    res = Parallel(n_jobs=args.jobs)(delayed(companion)(c, cl[c], args.n_sets, d[c]["tert"]) for c in con.use)
    res = dict(zip(con.use, res))
    log("companion sets done")
    comp = {}
    for kind in ("matched", "responsive"):
        comp[kind] = np.array([d_match(sum((res[c][kind][j] for c in con.use), []))[0] for j in range(args.n_sets)])
    t = dict(hypothesis="H185.1", D_match=Dm, p=p, null_mean=s["null_mean"], null_sd=s["null_sd"], null_q95=s["null_q95"], mde=s["mde"],
             n_perm=len(x), D_A175=D, n_bins=nb, n_bins_dropped=len(m.dropped), n_pairs_in_bins=npairs,
             dropped_bins=";".join(f"{c}:{k}({nl}/{nh})" for c, k, _, nl, nh in m.dropped))
    for kind, D_ in comp.items():
        ps = r75.pseudo_summary(D_, Dm)
        t |= {f"{kind}_q95": ps["q95"], f"{kind}_median": ps["median"], f"{kind}_percentile_of_D_match": ps["pct"],
              f"{kind}_exceeds_q95": ps["exceeds"], f"{kind}_n_sets": ps["n"]}
    t["clock_specific"] = bool(t["matched_exceeds_q95"] and t["responsive_exceeds_q95"])
    t["share_removed"] = (D - Dm) / D
    if p < 0.05 and Dm > 0:
        v = "passes: proliferating tumours follow their normal tissue less beyond having less clock signal"
        v += "; clock-specific" if t["clock_specific"] else "; not clock-specific"
    elif s["mde"] <= 0.1577:
        v = "fails, MDE <= 0.158: A175 pairing drop explained by less clock signal (swamped); wording: less clock signal, host following not shown to change"
    else:
        v = "not resolved"
    t["verdict"] = v
    pd.DataFrame([t]).to_csv(out / "tests.csv", index=False)
    pd.DataFrame({k + "_D_match": v for k, v in comp.items()}).to_csv(out / "companion_D_match.csv", index=False)

    # descriptives -------------------------------------------------------------------------------------------------
    desc = [dict(descriptive="D_A175", value=D), dict(descriptive="D_match", value=Dm),
            dict(descriptive="share of D removed (D - D_match)/D", value=t["share_removed"])]
    mn = Match(cs, d, cidx, con.use, amp_dict(cs, st, "zb"), 2)  # matching on normal amplitude
    nn = mn.null(args.n_perm)
    (Dn, nbn, npn), _ = mn.observed()
    xn = nn[np.isfinite(nn)]
    desc.append(dict(descriptive="D_match matched on normal-tissue amplitude |z_normal|", value=Dn,
                     p=float((1 + np.sum(xn >= Dn)) / (len(xn) + 1)), null_sd=float(xn.std(ddof=1)), mde=q.summarise(nn)["mde"],
                     n_bins=nbn, n_pairs=npn, share_removed=(D - Dn) / D, note=f"bins dropped {len(mn.dropped)}"))
    # amplitude distributions per arm (aggregates; amplitude relative to the cancer median over low+high pairs)
    for name, A in (("tumour |z|", amp_dict(cs, st, "za")), ("normal |z|", amp_dict(cs, st, "zb"))):
        for arm, al in ((0, "low G2M"), (2, "high G2M")):
            rel = np.concatenate([A[c][d[c]["tert"] == arm] / np.median(A[c][(d[c]["tert"] == 0) | (d[c]["tert"] == 2)]) for c in con.use])
            desc.append(dict(descriptive=f"{name} / cancer median, {al} arm (pool cancers)", value=float(np.median(rel)), n_pairs=len(rel),
                             q25=float(np.quantile(rel, .25)), q75=float(np.quantile(rel, .75)), mean=float(rel.mean())))
    # within-bin balance after tumour-amplitude matching: weighted mean of |za|/cancer median per arm
    for arm, al in ((0, "low"), (2, "high")):
        num = den = 0.0
        for c, k, idx, nl, nh in m.cells:
            med = np.median(m.amp[c][(m.lab[c] == 0) | (m.lab[c] == 2)])
            sel = idx[m.lab[c][idx] == arm]
            num += (m.amp[c][sel] / med).mean() * (nl + nh); den += nl + nh
        desc.append(dict(descriptive=f"tumour |z|/median, {al} arm, within matched bins (weights n_low+n_high)", value=num / den))
    pd.DataFrame(desc).to_csv(out / "descriptive.csv", index=False)
    print(pd.DataFrame([t]).T.to_string(), flush=True)
    print(pd.DataFrame(desc).round(4).to_string(), flush=True)
    log("tests stage done")


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["reproduce", "power", "tests"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-sets", type=int, default=200)
    args = ap.parse_args(argv)
    cli.setup(args)
    out = HERE / "output" / ("" if args.n_perm == 10_000 and args.n_sets == 200 else "debug")
    out.mkdir(parents=True, exist_ok=True)
    {"reproduce": stage_reproduce, "power": stage_power, "tests": stage_tests}[args.stage](args, out)


if __name__ == "__main__":
    main()
