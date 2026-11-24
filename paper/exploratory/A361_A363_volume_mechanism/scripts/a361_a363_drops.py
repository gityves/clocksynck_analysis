"""A361 / A363 gene-drop addendum: is the smaller core-loop coefficient carried by PER1, PER2 and NR1D1, or by the
precision lost when any three genes are dropped?

For every one of the 364 three-gene drops of the 14 genes: person volume on the remaining 11 genes (A360 rows), the GC
coefficient of A361's model and the inflammatory coefficient of A363's model. The trio drop is the core loop. Percentile
of the trio drop's coefficient among all drops and among the drops with similar singular values (A336 rule, within 20 %
of the core's); paired 14-minus-core difference of the inflammatory coefficient on 2,000 shared donor bootstraps (seed
[363, 4]); shrinkage predicted by the reliability loss (A336 part 2 half-host leverages, Spearman-Brown).

    python a361_a363_drops.py
Output: ../output/drops_coefs.csv, drops_tests.csv. Aggregates only.
"""
import itertools
import time

import numpy as np
import pandas as pd

import vm_common as V
from a361_a363_systemic import A360_MODEL, A363_MODEL, boot_coefs
from clocksync.config import limit_threads

TRIO = tuple(V.TRIO)
DROPS = list(itertools.combinations(V.WF.genes, 3))
LEV = {"14 genes": 0.765, "core": 0.714}   # A336 part 2 person leverages (half-host radii)


def sv(w):
    r = np.asarray(w.amp) / np.asarray(w.sd)
    return np.linalg.svd(np.c_[r * np.cos(w.psi), r * np.sin(w.psi)], compute_uv=False)


def main():
    V.quiet()
    limit_threads(4)
    t0 = time.time()
    pb, _, gc, ieg = V.a359.load()
    inf = V.inflammatory_genes(pb.columns)
    comps = sorted(pb.compartment.unique())
    per = V.persons(pb, comps, gc, ieg, inf)
    v14 = V.volume(pb, comps, V.WF)
    c14 = dict(gc=V.coefs(per, v14, A360_MODEL).gc, inf=V.coefs(per, v14, A363_MODEL).inf)
    rows, vols = [], {}
    for k, c in enumerate(DROPS):
        w = V.WF.drop(list(c))
        v = V.volume(pb, comps, w)
        if c == TRIO:
            vols["core"] = v
        s = sv(w)
        rows.append(dict(drop="+".join(c), trio=c == TRIO, sv1=s[0], sv2=s[1],
                         gc=V.coefs(per, v, A360_MODEL).gc, inf=V.coefs(per, v, A363_MODEL).inf))
        if k % 50 == 0:
            print(f"{k + 1}/{len(DROPS)} drops ({time.time() - t0:.0f}s)", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(V.OUT / "drops_coefs.csv", index=False)
    tr = d[d.trio].iloc[0]
    sim = (np.abs(d[["sv1", "sv2"]].to_numpy() / tr[["sv1", "sv2"]].to_numpy(float) - 1) <= 0.2).all(1)
    rel = lambda L: 2 * L / (1 + L)  # noqa: E731
    shrink = float(np.sqrt(rel(LEV["core"]) / rel(LEV["14 genes"])))
    bt, n = boot_coefs(per, {"14": v14, "core": vols["core"]}, A363_MODEL, [363, 4])
    diff = bt["14"][:, -1] - bt["core"][:, -1]
    out = []
    for term in ("gc", "inf"):
        # the effect is negative: the trio "carries" it if dropping the trio leaves a coefficient closer to 0 (higher)
        # than 90 % of the drops
        pct = float((d[term] < tr[term]).mean() * 100)
        pct_sim = float((d[term][sim] < tr[term]).mean() * 100) if sim.sum() >= 20 else np.nan
        r = dict(term=term, coef_14=c14[term], coef_core=float(tr[term]), predicted_core_from_precision=c14[term] * shrink,
                 drops_median=float(d[term].median()), drops_q10=float(d[term].quantile(0.10)),
                 drops_q90=float(d[term].quantile(0.90)), trio_pct_all=pct, n_similar=int(sim.sum()),
                 trio_pct_similar=pct_sim, carried_by_trio=bool(pct > 90))
        if term == "inf":
            r.update(diff_14_minus_core=c14[term] - float(tr[term]), diff_lo90=float(np.quantile(diff, 0.05)),
                     diff_hi90=float(np.quantile(diff, 0.95)), boot_persons=n)
        r["reading"] = ("carried by PER1, PER2, NR1D1" if r["carried_by_trio"] else
                        "not decided: precision or the removed genes (trio drop within the random drops)")
        out.append(r)
    t = pd.DataFrame(out)
    t["reliability_shrink_factor"] = shrink
    t.to_csv(V.OUT / "drops_tests.csv", index=False)
    print(t.round(4).T.to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
