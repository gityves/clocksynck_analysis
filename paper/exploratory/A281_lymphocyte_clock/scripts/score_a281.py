"""A281: do lymphocytes carry a clock of their own that the tissue does not share, or no readable clock?

Registered in paper/docs/registered_experiments.md (A281). Per compartment (>= 15 donors) of the healthy body atlases of
A272, A274 and A280: signed clock-gene coherence C across donors (`decompose.signed_coherence` on
log((count + 0.5) / depth) centred within stratum), one-sided p against per-gene permutations within strata. Per atlas
with a lymphocyte and a fibroblast-stromal compartment: both read on their shared donors, the fibroblast-stromal rows
thinned binomially to the lymphocyte compartment's median clock-gene count (20 replicates, 200 permutations each).

  H281.1  share of lymphocyte compartments with p < 0.05 >= 0.5
  H281.2  (read only if H281.1 fails) in atlases with thinned pass rate >= 0.8, lymphocytes fail in >= 2/3

    python score_a281.py [--design] [--jobs 6]
Outputs: paper/exploratory/A281_lymphocyte_clock/output/{design,compartments,pairs,tests}.csv
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.decompose import signed_coherence  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

spec = importlib.util.spec_from_file_location("a281_a280", PAPER / "exploratory/A280_celltype_replication2/scripts/score_a280.py")
A280 = importlib.util.module_from_spec(spec)
sys.modules["a281_a280"] = A280
spec.loader.exec_module(A280)  # A274 module with the anchor fallback; A150 / A272 through it
A274, A272, A150 = A280.A274, A280.A274.A272, A280.A274.A150

AMENDMENT = "A281_lymphocyte_clock"
SEED = 281
ATLASES = {"lung_hlca": "lung_hlca_genome", "heart": "heart_global_genome", "organ_breast": "organ_breast_genome",
           "organ_colon": "organ_colon_genome", "organ_muscle": "organ_muscle_genome",
           "organ_optic_nerve": "organ_optic_nerve_genome", "organ_eye_anterior": "organ_eye_anterior_genome",
           "kidney_sc": "kidney_sc_genome", "kong_ileum": "kong_ileum_genome", "kong_colon": "kong_colon_genome",
           "ibd_ileum": "ibd_ileum_genome", "ibd_rectum": "ibd_rectum_genome"}
MIN_DONORS, N_PERM, N_THIN, N_PERM_THIN, N_PSEUDO = 15, 1000, 20, 200, 200
WF = cs.gtex()
CLOCK = list(WF.genes)


def y_of(t, genes):
    """log((count + 0.5) / depth), centred within stratum."""
    y = np.log((t[genes].to_numpy(float) + 0.5) / t.depth.to_numpy(float)[:, None])
    for idx in t.groupby("stratum").indices.values():
        y[idx] -= y[idx].mean(0)
    return y


def perm_test(y, groups, n_perm, rng, wf=WF):
    obs = signed_coherence(y, wf)
    null = np.empty(n_perm)
    for k in range(n_perm):
        z = y.copy()
        for i in groups:
            for j in range(y.shape[1]):
                z[i, j] = y[rng.permutation(i), j]
        null[k] = signed_coherence(z, wf)
    return obs, float((1 + np.sum(null >= obs)) / (1 + n_perm))


def groups_of(t):
    return [np.asarray(i) for i in t.reset_index(drop=True).groupby("stratum").indices.values()]


def thin(t, target, rng):
    """Rows thinned binomially so that their clock-gene sum is about `target` (rows below it are kept)."""
    t = t.copy()
    s = t[CLOCK].sum(1).to_numpy(float)
    p = np.minimum(1.0, target / np.maximum(s, 1))
    t[CLOCK] = rng.binomial(t[CLOCK].to_numpy(np.int64), p[:, None])
    t["depth"] = t.depth.to_numpy(float) * p
    return t


def compartment_task(n, comp, t, seed):
    warnings.filterwarnings("ignore")
    t = t.reset_index(drop=True)
    c, p = perm_test(y_of(t, CLOCK), groups_of(t), N_PERM, np.random.default_rng(seed))
    return dict(atlas=n, compartment=comp, family=A150.FAM.family(comp), donors=t.donor.nunique(),
                clock_count_median=float(t[CLOCK].sum(1).median()), C=c, p=p)


def pair_task(n, ly, sv, tl, ts, pseudo, seed):
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(seed)
    d = sorted(set(tl.donor) & set(ts.donor))
    tl = tl[tl.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    ts = ts[ts.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    target = float(tl[CLOCK].sum(1).median())
    cl, pl = perm_test(y_of(tl, CLOCK), groups_of(tl), N_PERM, rng)
    cs_, ps_ = perm_test(y_of(ts, CLOCK), groups_of(ts), N_PERM, rng)
    thin_p, thin_c = [], []
    for _ in range(N_THIN):
        tt = thin(ts, target, rng)
        c, p = perm_test(y_of(tt, CLOCK), groups_of(tt), N_PERM_THIN, rng)
        thin_c.append(c), thin_p.append(p)
    pc = np.array([signed_coherence(y_of(tl, list(w.genes)), w) for w in pseudo])
    return dict(atlas=n, lymphocyte=ly, stromal=sv, shared_donors=len(d), lymph_clock_count_median=target,
                stromal_clock_count_median=float(ts[CLOCK].sum(1).median()), C_lymph=cl, p_lymph=pl, C_stromal=cs_,
                p_stromal=ps_, C_stromal_thinned_median=float(np.median(thin_c)),
                thinned_pass_rate=float(np.mean(np.array(thin_p) < 0.05)),
                lymph_pct_among_matched=float(np.mean(pc[np.isfinite(pc)] < cl)), matched_C_median=float(np.nanmedian(pc)))


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--design", action="store_true")
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(AMENDMENT)
    out.mkdir(parents=True, exist_ok=True)
    A150.ATLASES.update(ATLASES)
    tabs, design, pairs_in = {}, [], []
    for k, n in enumerate(ATLASES):
        A = A150.load_atlas(n, None, 0)
        pb = A["pb"].reset_index(drop=True)
        nd = pb.groupby("compartment").donor.nunique()
        comps = sorted(nd.index[nd >= MIN_DONORS])
        fam = {c: A150.FAM.family(c) for c in comps}
        ly = [c for c in comps if fam[c] == "lymphocyte"]
        fb = [c for c in comps if fam[c] == "fibroblast_stromal"]
        tabs[n] = (A, pb, comps)
        row = dict(atlas=n, donors=pb.donor.nunique(), compartments=len(comps), lymphocyte=";".join(ly),
                   fibroblast_stromal=";".join(fb))
        if ly and fb:
            l, f = max(ly, key=lambda c: nd[c]), max(fb, key=lambda c: nd[c])
            sh = len(set(pb.donor[pb.compartment == l]) & set(pb.donor[pb.compartment == f]))
            row.update(pair_lymphocyte=l, pair_stromal=f, shared_donors=sh)
            if sh >= MIN_DONORS:
                pairs_in.append((k, n, l, f))
        design.append(row)
    design = pd.DataFrame(design)
    design.to_csv(out / "design.csv", index=False)
    print(design.to_string())
    if args.design:
        return

    tasks = [delayed(compartment_task)(n, c, pb[pb.compartment == c], [SEED, k, j])
             for k, (n, (A, pb, comps)) in enumerate(tabs.items()) for j, c in enumerate(comps)]
    comp = pd.DataFrame(Parallel(n_jobs=args.jobs)(tasks))
    comp.to_csv(out / "compartments.csv", index=False)
    ptasks = []
    for k, n, l, f in pairs_in:
        A, pb, _ = tabs[n]
        pseudo = A272.pseudo_sets(A, k, N_PSEUDO, "matched")
        cols = sorted({g for w in pseudo for g in w.genes} | set(CLOCK))
        sub = pb[["donor", "compartment", "stratum", "n_cells", "depth"] + cols]
        ptasks.append(delayed(pair_task)(n, l, f, sub[sub.compartment == l], sub[sub.compartment == f], pseudo,
                                         [SEED, 100 + k]))
    pairs = pd.DataFrame(Parallel(n_jobs=args.jobs)(ptasks))
    pairs.to_csv(out / "pairs.csv", index=False)
    share = float(np.mean(pairs.p_lymph < 0.05))
    powered = pairs[pairs.thinned_pass_rate >= 0.8]
    h2_share = float(np.mean(powered.p_lymph >= 0.05)) if len(powered) else np.nan
    tests = pd.DataFrame([dict(atlases=len(pairs), lymph_pass=int((pairs.p_lymph < 0.05).sum()), H281_1_share=share,
                               H281_1_pass=share >= 0.5, powered_atlases=len(powered),
                               powered_lymph_fail=int((powered.p_lymph >= 0.05).sum()) if len(powered) else 0,
                               H281_2_share_fail=h2_share,
                               H281_2_pass=(bool(h2_share >= 2 / 3) if share < 0.5 and len(powered) else None))])
    tests.to_csv(out / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(pairs.round(3).to_string())
    print(comp.groupby("family").agg(n=("C", "size"), C_median=("C", "median"), share_p05=("p", lambda p: float(np.mean(p < 0.05)))).round(3))
    print(tests.T.to_string())


if __name__ == "__main__":
    main()
