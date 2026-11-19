"""A333 (registered, check N1): without behaviour, do a tissue's daily outputs follow its clock?

Per A109 body tissue: fast deaths (Hardy 1-2) split once into halves A (selection) and B (reference); V = the ventilated
post-mortem donors with a living brain (Hardy 0, cohort Postmortem; 79). Output program = non-clock, non-responsive genes
(median CPM >= 5) called clock-following in half A by the A91 model (internal p < 0.01 given wall time, wall-time
p > 0.05; D2 covariates; host from the other body tissues). Per group g in {B, V}, on the group's own leave-one-out host:
L_out = median over program genes of |coefficient of log CPM on the host state| (D2 covariates); L_clock = |z| of the
tissue's clock-plane state on the same host. Phi = (L_out,V / L_out,B) / (L_clock,V / L_clock,B).

    python score_a333.py --stage power   # null construction only: V labels permuted among B and V donors (200)
    python score_a333.py --stage score   # observed Phi, bootstrap, companions, matching (only if the gate passed)
Outputs (`../output/`): program.csv, power.csv, power_perm.csv | (score stage, not yet written). Aggregates only
(GTEx controlled access: per-donor values stay in memory).
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_3" / "scripts" / "gtex_cues"))
sys.path.insert(0, str(PAPER / "exploratory" / "A132_brain_death" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE, RESPONSIVE_GENES  # noqa: E402
from clocksync._core.responsive_broad import RESPONSIVE_BROAD  # noqa: E402
from clocksync.cells.dispersion import _projector  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

import _shared as S  # noqa: E402
import compute_brain_death as B  # noqa: E402

OUT = HERE.parents[1] / "output"
N_PERM, MIN_FAST, MIN_V, MIN_PROGRAM = 200, 40, 30, 20
EXCLUDE = set(CLOCK_LIKE) | set(RESPONSIVE_GENES) | set(RESPONSIVE_BROAD)
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:6.0f}s] {m}", flush=True)


def groups():
    """Donor sets: fast halves A / B and V (living-brain ventilated post-mortem)."""
    ph = gtex.phenotypes(["COHORT", "DTHHRDY"])
    pb, _ = S.clock_table(hardy=(0, 1, 2))
    donors = sorted(pb.donor.unique())
    p = ph.reindex(donors)
    fast = np.array([d for d in donors if p.DTHHRDY.get(d) in (1, 2)])
    V = {d for d in donors if p.DTHHRDY.get(d) == 0 and p.COHORT.get(d) == B.PM}
    fast = fast[np.random.default_rng([333, 1]).permutation(len(fast))]
    return set(fast[: len(fast) // 2]), set(fast[len(fast) // 2:]), V, pb


def complex_resid(M, Z):
    D = np.c_[np.ones(len(M)), Z]
    return M - D @ np.linalg.lstsq(D, M, rcond=None)[0]


class Tissue:
    """Rows of one tissue (one sample per donor), program genes, clock-plane projector."""

    def __init__(self, tis, A, Bset, V, pb):
        self.tis = tis
        self.others = [t for t in S.BODY if t != tis]
        Y, genes, ob = S.expression(tis, hardy=(0, 1, 2), tod=False, median_cpm=None)
        self.subj = ob.subject.to_numpy()
        self.ob = ob.reset_index(drop=True)
        wf = cs.gtex()
        gi = {g: i for i, g in enumerate(genes)}
        present = [g for g in wf.genes if g in gi]
        self.ck = [gi[g] for g in present]
        self.L = _projector(wf.subset(present).X, np.ones(len(present)))
        inA = np.isin(self.subj, list(A))
        self.n_fast_B, self.n_V = int(np.isin(self.subj, list(Bset)).sum()), int(np.isin(self.subj, list(V)).sum())
        # program: selection in half A (rows with a time of death), A91 model
        a = inA & ob.tod.notna().to_numpy()
        expr = np.median(np.exp(Y[a]) - 0.5, 0) >= 5
        cand = expr & ~np.isin([str(g) for g in genes], list(EXCLUDE))
        hA = S.host_states(pb, self.others, donors=A)
        m = a & np.isin(self.subj, hA.index)
        f = hA[["x", "y"]].loc[self.subj[m]].to_numpy()
        f = f - f.mean(0)
        t = ob.tod.to_numpy(float)[m] % 24
        E = np.c_[np.cos(S.W * t), np.sin(S.W * t)]
        Zc = S.d2_covariates(ob[m])
        Zc = Zc[:, Zc.std(0) > 0]
        D = np.c_[np.ones(m.sum()), Zc, f, E]
        k = D.shape[1]
        Yc = Y[m][:, cand]
        _, p_int = S.partial_f(Yc, D, [k - 4, k - 3])
        _, p_ext = S.partial_f(Yc, D, [k - 2, k - 1])
        sel = (p_int < 0.01) & (p_ext > 0.05)
        self.prog = np.flatnonzero(cand)[sel]
        self.n_selA = int(m.sum())
        self.Yp = Y[:, self.prog]
        self.Yk = Y[:, self.ck]

    def loadings(self, donors, host):
        """L_out and L_clock of the rows of `donors` on `host` (x, y by donor)."""
        m = np.isin(self.subj, list(donors)) & np.isin(self.subj, host.index)
        h = host[["x", "y"]].loc[self.subj[m]].to_numpy()
        Zc = S.d2_covariates(self.ob[m])
        Zc = Zc[:, Zc.std(0) > 0]
        hr = complex_resid(h, Zc)
        hc = hr[:, 0] + 1j * hr[:, 1]
        yp = complex_resid(self.Yp[m], Zc)
        b = np.linalg.lstsq(hr, yp, rcond=None)[0]                 # 2 x genes
        L_out = float(np.median(np.hypot(b[0], b[1])))
        s = complex_resid(self.Yk[m], Zc) @ self.L.T               # tissue clock-plane state
        sc = s[:, 0] + 1j * s[:, 1]
        z = np.sum(sc * np.conj(hc)) / np.sum(np.abs(hc) ** 2)
        return L_out, float(abs(z)), int(m.sum())


def log_phi(T, Bd, Vd, pb):
    out = []
    for d in (Bd, Vd):
        host = S.host_states(pb, T.others, donors=d)
        out.append(T.loadings(d, host))
    (ob_, cb, _), (ov, cv, _) = out
    return float(np.log((ov / ob_) / (cv / cb)))


def power(jobs):
    A, Bd, V, pb = groups()
    log(f"donors: half A {len(A)}, half B {len(Bd)}, V {len(V)}")
    tis, rows = [], []
    for t in S.BODY:
        T = Tissue(t, A, Bd, V, pb)
        ok = T.n_fast_B + len(A & set(T.subj)) >= MIN_FAST and T.n_V >= MIN_V and len(T.prog) >= MIN_PROGRAM
        rows.append(dict(tissue=t, fast_B=T.n_fast_B, V=T.n_V, selection_rows=T.n_selA, program_genes=len(T.prog),
                         estimable=bool(ok)))
        log(str(rows[-1]))
        if ok:
            tis.append(T)
    pd.DataFrame(rows).to_csv(OUT / "program.csv", index=False)
    pool = np.array(sorted(Bd | V))
    nV = len(V)

    def perm(s):
        r = np.random.default_rng([333, 2, s])
        lab = r.permutation(len(pool))
        Vp, Bp = set(pool[lab[:nV]]), set(pool[lab[nV:]])
        return [dict(perm=s, tissue=T.tis, log_phi=log_phi(T, Bp, Vp, pb)) for T in tis]
    res = Parallel(n_jobs=jobs)(delayed(perm)(s) for s in range(N_PERM))
    pp = pd.DataFrame([x for r in res for x in r])
    pp.to_csv(OUT / "power_perm.csv", index=False)
    pw = pp.groupby("tissue").log_phi.agg(null_mean="mean", se="std").reset_index()
    pw["MDE"] = 2.49 * pw.se
    pw["MDE_ratio"] = np.exp(-pw.MDE)
    med = float(pw.MDE.median())
    pw.loc[len(pw)] = dict(tissue="median over tissues", MDE=med, MDE_ratio=float(np.exp(-med)))
    pw["gate"] = None
    pw.loc[pw.tissue == "median over tissues", "gate"] = bool(med <= np.log(1 / 0.7))
    pw.to_csv(OUT / "power.csv", index=False)
    log(pw.round(3).to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "score":
        pw = pd.read_csv(OUT / "power.csv")
        if str(pw.loc[pw.tissue == "median over tissues", "gate"].iloc[0]) != "True":
            raise SystemExit("power gate not passed: not decidable at this design; not scored")
        raise SystemExit("score stage not written yet: write it (as registered) before scoring")
    power(a.jobs)
    log("done")


if __name__ == "__main__":
    main()
