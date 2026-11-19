"""B1 power scout (prompt cell_route_2026-11-19, Task 4): does the GTEx fast-death organ arrangement agree with the
baboon's (A36 timed tissue phases) on the shared A109 tissues? Fisher-Lee rho over the shared pairs, tissue-label
permutation null (exact over 6! orders). Groups' rho reported too (descriptive). Aggregates only.
Output: ../output/b1_scout.csv"""
import itertools
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a317 as A  # noqa: E402
from score_a317 import B, H, ST, gtex  # noqa: E402

bab = pd.read_csv(HERE.parents[3] / "results/figure_1/output/species/a36_pairs.csv").set_index("gtex")
shared = [t for t in H.BODY if t in bab.index]
T = H.Tissues(200, 0)
T.build_states(np.random.default_rng(1), 0, real=True)
ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
coh = ph.COHORT.to_numpy()
fast, h4, h0 = (T.group_donors(g) for g in ("fast", "H4", "vent"))
od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
draws = [np.random.default_rng(A.SEED).choice(od, len(pm), replace=False)]
Z, M = A.raw(T, "clock", 0, T.R)
idx = [H.BODY.index(t) for t in shared]
pairs = list(itertools.combinations(range(len(idx)), 2))
mask = np.array([(a in idx and b in idx) for a, b in H.PAIRS])
# baboon pair phase: phase(a) - phase(b) in the same sign as arg E[z_a conj z_b]; the sign of the human/baboon convention
# is checked by reading both orientations
phi = bab.loc[shared, "baboon_g_phi_h"].to_numpy() / H.HOURS
bab_pairs = np.array([phi[i] - phi[j] for i, j in pairs])
rows = []
groups = {"fast": fast, "hardy4": h4, "postmortem_vent": pm, "organ_donors_all": od}
for g, d in groups.items():
    th = B.phase(Z, M, d, mask)
    for sign in (1, -1):
        r = ST.fisher_lee(th, sign * bab_pairs)
        null = [ST.fisher_lee(th, sign * np.array([phi[p[i]] - phi[p[j]] for i, j in pairs]))
                for p in itertools.permutations(range(len(idx)))]
        rows.append(dict(group=g, sign=sign, tissues=len(idx), pairs=len(pairs), n_donors=len(d), rho=float(r),
                         p_perm=float(np.mean(np.array(null) >= r))))
out = pd.DataFrame(rows)
out.to_csv(A.OUT / "b1_scout.csv", index=False)
print(shared)
print(out.round(3).to_string(index=False))
