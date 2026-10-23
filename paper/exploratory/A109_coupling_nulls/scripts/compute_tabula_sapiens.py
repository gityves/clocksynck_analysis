"""A109b (secondary, not in the confirmatory family): replication in Tabula Sapiens stromal-vascular pseudobulks
(brain-dead ventilated organ donors; 12 donors, 9 organs with >= 5 donors, 49 donor x organ pseudobulks).

Registered in paper/docs/registered_experiments.md (A109b and its note). Two statistics:
  C against the pseudo-clock null (as H109.1). C = mean over organ pairs with >= `min_shared` shared donors of the
      Poisson-corrected |complex correlation| (A92 form, centred over the shared donors), the modulus clipped at 1
      (with 3-5 donors the corrected denominators can fall below the covariance); null: 200 expression-matched
      14-gene sets (`clocksync.paired.pseudo_waveforms`, matched per organ on mean log CPM, clock, clock-like and D6
      clock-output genes excluded). Statistic C minus null median, one-sided add-one p. Run at min_shared 3 (primary,
      all 36 pairs) and 5 (7 pairs).
  Offset preservation: phase arg E[z_a conj z_b] of the 6 pairs of the 4 organs shared with GTEx (adipose ->
      Adipose - Subcutaneous, heart -> Heart - Atrial Appendage, musculature -> Muscle - Skeletal, respiratory system
      -> Lung) against the GTEx fast-death pair offsets (A109 pair_relations.csv); Fisher-Lee correlation, tissue
      labels permuted jointly (4 organs: 24 permutations, exact one-sided p, minimum 1/24 = 0.042).
--mde computes the pseudo-clock null of C and the power of the offset test only; it never computes C on the real
clock genes. Without --mde the scoring is run (after the MDE is committed).

    python compute_tabula_sapiens.py [--mde] [--subsample N] [--seed S]
Outputs (paper/exploratory/A109_coupling_nulls/output/, .../subsample/): a109b_c_null.csv, a109b_offsets.csv,
a109b_power.csv (--mde).
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from common import cli  # noqa: E402
from common.paths import PSEUDOBULK, exploratory_dir  # noqa: E402

FOLDER = "A109_coupling_nulls"
P, WF = H.P, H.WF
GTEX_OF = {"adipose tissue": "Adipose - Subcutaneous", "heart": "Heart - Atrial Appendage",
           "musculature": "Muscle - Skeletal", "respiratory system": "Lung"}


class Organs:
    def __init__(self, n_sets: int, seed: int, subsample: int | None):
        pb = pd.read_parquet(PSEUDOBULK / "tabula_sapiens_organs.parquet")
        if subsample:
            keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
            pb = pb[pb.donor.isin(keep)]
        pb = pb.assign(organ=pb.compartment.str.split(":").str[0])
        self.organs = sorted(pb.organ.unique())
        self.donors = sorted(pb.donor.unique())
        gi = {d: i for i, d in enumerate(self.donors)}
        clock = list(WF.genes)
        self.R, self.Z, self.V, self.sets = {}, {"clock": {}, "pseudo": {}}, {"clock": {}, "pseudo": {}}, {}
        for k, o in enumerate(self.organs):
            x = pb[pb.organ == o].reset_index(drop=True)
            dep = x.depth.to_numpy(float)
            r = np.full(len(self.donors), -1)
            r[[gi[d] for d in x.donor]] = np.arange(len(x))
            self.R[o] = r
            genes = x.drop(columns=["donor", "compartment", "stratum", "n_cells", "depth", "organ"])
            lm = pd.Series(np.log(genes.to_numpy(float) / dep[:, None] * 1e6 + 0.5).mean(0), index=genes.columns)
            lm = lm[~lm.index.isin(H.CLOCK_OUTPUT - set(clock))]
            lm = lm[~lm.index.duplicated()]
            wfs = P.pseudo_waveforms(lm, n_sets, np.random.default_rng(1000 * (seed + 1) + k), WF)
            self.Z["clock"][o], self.V["clock"][o] = H._state_set(x[clock], [WF], dep)
            z, v = [], []
            for w in wfs:
                xx = genes[list(w.genes)].astype(float)
                b, vv = P.states(xx, P.projector(xx, w, counts=True), w, depth=dep)
                z.append(b[:, 0] + 1j * b[:, 1])
                v.append(vv)
            self.Z["pseudo"][o], self.V["pseudo"][o] = np.array(z), np.array(v)

    def pairs(self, min_shared):
        return [(a, b) for a, b in itertools.combinations(self.organs, 2)
                if ((self.R[a] >= 0) & (self.R[b] >= 0)).sum() >= min_shared]

    def corr(self, kind, min_shared):
        """(K, n_pairs) clipped |ccorr| of the pairs with >= min_shared shared donors."""
        out = []
        for a, b in self.pairs(min_shared):
            m = (self.R[a] >= 0) & (self.R[b] >= 0)
            ia, ib = self.R[a][m], self.R[b][m]
            A, B = self.Z[kind][a][:, ia], self.Z[kind][b][:, ib]
            A, B = A - A.mean(1, keepdims=True), B - B.mean(1, keepdims=True)
            k = m.sum() - 1
            ea = (abs(A) ** 2).sum(1) / k - self.V[kind][a][:, ia].mean(1)
            eb = (abs(B) ** 2).sum(1) / k - self.V[kind][b][:, ib].mean(1)
            c = abs((A * B.conj()).sum(1) / k) / np.sqrt(np.maximum(ea, 1e-12) * np.maximum(eb, 1e-12))
            out.append(np.minimum(c, 1.0))
        return np.array(out).T

    def C(self, kind, min_shared):
        c = self.corr(kind, min_shared)
        return np.nanmean(c, axis=1) if c.ndim == 2 else np.full(len(self.Z[kind][self.organs[0]]), np.nan)


def pval(obs, null):
    null = np.asarray(null)
    return float((1 + np.sum(null >= obs)) / (1 + len(null)))


def gtex_reference():
    """GTEx fast-death offsets (h, A minus B) of the 6 pairs of the 4 shared organs, in the order of `pairs4`."""
    rel = pd.read_csv(exploratory_dir(FOLDER) / "pair_relations.csv")
    rel = rel[rel.group == "fast"].set_index(["tissue_a", "tissue_b"]).offset_h
    organs = sorted(GTEX_OF)
    ref, flip = [], []
    for a, b in itertools.combinations(organs, 2):
        ga, gb = GTEX_OF[a], GTEX_OF[b]
        ref.append(rel[(ga, gb)] if (ga, gb) in rel.index else -rel[(gb, ga)])
    return organs, np.array(ref) / H.HOURS


def offset_matrix(O: Organs, organs):
    """4 x 4 complex E[z_a conj z_b] over shared donors (centred over the shared donors), NaN with < 3 donors."""
    Hm = np.full((len(organs),) * 2, np.nan + 0j)
    for i, a in enumerate(organs):
        for j, b in enumerate(organs):
            if i == j:
                continue
            m = (O.R[a] >= 0) & (O.R[b] >= 0)
            if m.sum() < 3:
                continue
            A, B = O.Z["clock"][a][0][O.R[a][m]], O.Z["clock"][b][0][O.R[b][m]]
            Hm[i, j] = np.sum((A - A.mean()) * np.conj(B - B.mean())) / (m.sum() - 1)
    return Hm


def fl_test(theta_ref, Hm):
    """Fisher-Lee of the 6 pair phases against the reference; the 24 joint permutations of the 4 organ labels."""
    iu = np.triu_indices(4, 1)
    if np.isnan(Hm[iu]).any():  # a pair with fewer than 3 shared donors: not estimable
        return np.nan, np.nan, np.array([np.nan])
    obs = ST.fisher_lee(theta_ref, np.angle(Hm)[iu])
    null = np.array([ST.fisher_lee(theta_ref, np.angle(Hm[np.ix_(p, p)])[iu]) for p in itertools.permutations(range(4))])
    return obs, float(np.sum(null >= obs - 1e-12) / 24), null


def power(O: Organs, organs, theta_ref, rng, n_sim):
    """Power of the offset test: rank-1 clocks with GTEx-fast offsets (the leading eigenvector of the reference phases
    is rebuilt from the 6 pair phases), shared donor pattern of the data, coupling c; exact permutation p."""
    iu = np.triu_indices(4, 1)
    th = np.zeros((4, 4))
    th[iu] = theta_ref
    th = th - th.T
    delta = np.angle(np.exp(1j * th).sum(1))  # offset of each organ from the mean of the others
    rows = []
    for c in (0.3, 0.45, 0.6, 0.8, 0.95):
        rs, ps = [], []
        for _ in range(n_sim):
            phi = {d: rng.uniform(0, 2 * np.pi) for d in O.donors}
            Z = {}
            for i, o in enumerate(organs):
                z = np.zeros(len(O.donors), complex)
                for j, d in enumerate(O.donors):
                    if O.R[o][j] >= 0:
                        xi = (rng.standard_normal() + 1j * rng.standard_normal()) / np.sqrt(2)
                        z[j] = np.sqrt(c) * np.exp(1j * (phi[d] + delta[i])) + np.sqrt(1 - c) * xi
                Z[o] = z
            Hm = np.full((4, 4), np.nan + 0j)
            for i, a in enumerate(organs):
                for j, b in enumerate(organs):
                    m = (O.R[a] >= 0) & (O.R[b] >= 0)
                    if i != j and m.sum() >= 3:
                        A, B = Z[a][m], Z[b][m]
                        Hm[i, j] = np.sum((A - A.mean()) * np.conj(B - B.mean()))
            r, p, _ = fl_test(theta_ref, Hm)
            rs.append(r)
            ps.append(p)
        ps = np.array(ps)
        rows.append(dict(coupling=c, mean_rho=float(np.nanmean(rs)), power_p05=float(np.mean(ps < 0.05)),
                         power_p0417=float(np.mean(ps <= 1 / 24 + 1e-9)), min_p=1 / 24, n_sim=n_sim))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--mde", action="store_true", help="null of C and power only (no real C)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    out = exploratory_dir(FOLDER, quick)
    rng = np.random.default_rng(1109 + args.seed)
    n_sets = 20 if quick else 200
    O = Organs(n_sets, args.seed, args.subsample)
    print(len(O.donors), "donors", len(O.organs), "organs; pairs >= 3 / >= 5 shared:", len(O.pairs(3)), len(O.pairs(5)),
          flush=True)
    rows = []
    for m in (3, 5):
        if not O.pairs(m):
            continue
        nul = O.C("pseudo", m)
        row = dict(min_shared=m, null_median=np.median(nul), null_sd=nul.std(ddof=1), null_q95=np.quantile(nul, 0.95),
                   mde_p05=np.quantile(nul, 0.95) - np.median(nul), n_sets=len(nul), n_pairs=len(O.pairs(m)),
                   n_donors=len(O.donors))
        if not args.mde:
            c = float(O.C("clock", m)[0])
            row.update(C_clock=c, excess=c - np.median(nul), p_one_sided=pval(c, nul))
        rows.append(row)
    cn = pd.DataFrame(rows)
    cn.to_csv(out / "a109b_c_null.csv", index=False)
    print(cn.round(4).to_string(index=False), flush=True)
    organs, theta_ref = gtex_reference()
    organs = [o for o in organs if o in O.organs]
    if len(organs) < 4:
        print("fewer than 4 shared organs in the subsample: offset test skipped")
        return
    if args.mde:
        pw = power(O, organs, theta_ref, rng, 30 if quick else 500)
        pw.to_csv(out / "a109b_power.csv", index=False)
        print(pw.round(3).to_string(index=False))
        return
    Hm = offset_matrix(O, organs)
    r, p, null = fl_test(theta_ref, Hm)
    iu = np.triu_indices(4, 1)
    pr = [dict(kind="pair", organ_a=organs[a], organ_b=organs[b], n_shared=int(((O.R[organs[a]] >= 0) &
              (O.R[organs[b]] >= 0)).sum()), offset_h=float(np.angle(Hm[a, b]) * H.HOURS),
              gtex_fast_offset_h=float(theta_ref[k] * H.HOURS)) for k, (a, b) in enumerate(zip(*iu))]
    pr.append(dict(kind="summary", rho_fisher_lee=r, p=p, n_pairs=6, n_perm=24, min_p=1 / 24, null_median=np.median(null)))
    pd.DataFrame(pr).to_csv(out / "a109b_offsets.csv", index=False)
    print(pd.DataFrame(pr).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
