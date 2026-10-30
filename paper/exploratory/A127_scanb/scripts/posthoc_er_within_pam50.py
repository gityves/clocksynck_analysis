"""A127 POST HOC (2026-10-30, chosen after H127.1-2 were scored; handoff item 2): is the ER- coherence loss present within
PAM50 subtype? ER is nearly collinear with PAM50 (LumA 9 / 1,633 ER-, LumB 1 / 727), so only Basal and Her2 carry both arms.

T_pam = mean over {Basal, Her2} of coherence(ER-) - coherence(ER+), ER labels permuted within subtype (one-sided, loss).
T_pam_prolif = the same within subtype x proliferation tertile (cohort tertiles), cells with >= 10 tumours per ER arm.
Specificity: the 200 variance-matched pseudo sets of run_a127.py (`pseudo_sets_variance_matched.csv`), observed labels.
Descriptive: per-subtype differences, proliferation score by ER within subtype.
    python -I posthoc_er_within_pam50.py [--nperm 2000]
Writes ../output/posthoc_er_within_pam50.csv, posthoc_er_within_pam50_cells.csv.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a127_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--nperm", type=int, default=2000)
args = ap.parse_args()
LN2, SUB, MIN_ARM = np.log(2), ("Basal", "Her2"), 10
rng = np.random.default_rng(1271)
m, X = L.load()
er, pam = m.er.to_numpy(), m.pam50.to_numpy()
prol_t = L.tertile(m.prolif.to_numpy())
Yc = LN2 * X[L.GENES].to_numpy(float)
sets = pd.read_csv(L.OUT / "pseudo_sets_variance_matched.csv").to_numpy()


def groups(by_prolif: bool) -> list:
    out = []
    for s in SUB:
        for t in (range(3) if by_prolif else [None]):
            i = np.flatnonzero((pam == s) & ((prol_t == t) if by_prolif else True))
            if min((er[i] == 0).sum(), (er[i] == 1).sum()) >= MIN_ARM:
                out.append(i)
    return out


def T(Y, wf, G, perms=None):
    d = []
    for k, i in enumerate(G):
        lab = er[i] if perms is None else er[i][perms[k]]
        d.append(L.coh(Y, i[lab == 0], wf) - L.coh(Y, i[lab == 1], wf))
    return float(np.mean(d)), d


rows, cells = [], []
for name, G in (("ER- minus ER+ within PAM50 (Basal, Her2)", groups(False)),
                ("ER- minus ER+ within PAM50 x proliferation tertile", groups(True))):
    obs, per = T(Yc, L.WF, G)
    null = np.array([T(Yc, L.WF, G, [rng.permutation(len(i)) for i in G])[0] for _ in range(args.nperm)])
    tp = np.array([T(LN2 * X[list(s)].to_numpy(float), L.wf_for(s), G)[0] for s in sets])
    rows.append(dict(analysis=name, T=obs, p_one_lower=cs.stats.perm_p(null, obs, "less", tie=0.0), null_sd=null.std(ddof=1),
                     MDE_80pct=(norm.ppf(0.95) + norm.ppf(0.80)) * null.std(ddof=1), n_cells=len(G),
                     n_tumours=int(sum(len(i) for i in G)), pseudo_percentile=100 * np.mean(tp <= obs),
                     pseudo_q05=np.quantile(tp, 0.05), n_perm=args.nperm, n_pseudo=len(tp)))
    for i, dv in zip(G, per):
        cells.append(dict(analysis=name, pam50=pam[i[0]], n_er_neg=int((er[i] == 0).sum()), n_er_pos=int((er[i] == 1).sum()),
                          diff=dv, prolif_er_neg=m.prolif.to_numpy()[i][er[i] == 0].mean(),
                          prolif_er_pos=m.prolif.to_numpy()[i][er[i] == 1].mean()))
res = pd.DataFrame(rows)
res.to_csv(L.OUT / "posthoc_er_within_pam50.csv", index=False)
pd.DataFrame(cells).to_csv(L.OUT / "posthoc_er_within_pam50_cells.csv", index=False)
print(res.round(4).to_string())
print(pd.DataFrame(cells).round(3).to_string())
