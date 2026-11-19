"""A335 part 1c: A331 (person loudness, MuTHER fat x skin) rerun with the core loop.

Reuses `score_a331.py` (cohort loader, leverage, permutation null with the same rng). Companions: A331's matched (tag 3)
and responsive (tag 4) sets, 200 each, drawn from the cohort waveform and restricted to the core positions (A329 rule).
The power gate is inherited (not recomputed).
    python a335_1c_muther.py [--check]     # --check: 3 sets, 20 permutations, assertions only
Output: ../output/a335_1c_tests.csv, a335_1c_companions.csv.
"""
import argparse
import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_2" / "scripts" / "backbone"))
spec = importlib.util.spec_from_file_location("a331", PAPER / "exploratory/P1_person_loudness/scripts/score_a331.py")
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)
from clocksync import elasticity as E  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
DROP = ["PER1", "PER2", "NR1D1"]


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    m = M.a255()
    c, persons = m.load_muther(False)
    both = persons[(persons.has_fat == "True") & (persons.has_skin == "True")].person
    idx = np.array([i in set(both) for i in c.ids])
    wf = c.wf
    core = wf.drop(DROP)
    core_idx = [k for k, g in enumerate(wf.genes) if g in core.genes]
    assert len(core_idx) == len(core.genes) == len(wf.genes) - 3
    n_sets = 3 if chk else m.N_SETS
    lev = lambda w: E.radius_covariance(*[c.Z(w)[t][idx] for t in ("fat", "skin")])["leverage"]  # noqa: E731
    pw = pd.read_csv(OUT.parents[1] / "P1_person_loudness/output/power.csv", index_col=0)["value"]
    assert str(pw["gate"]) == "True"
    rows = []
    for kind, tag in (("matched", 3), ("responsive", 4)):
        for j, w in enumerate(c.pseudo(kind, n_sets, tag)):
            wc = Waveform(tuple(w.genes[k] for k in core_idx), core.amp, core.psi, None, None, f"core{j}")
            Zf, Zc = c.Z(w), c.Z(wc)
            rows.append(dict(kind=kind, set=j, leverage_14gene=E.radius_covariance(Zf["fat"][idx], Zf["skin"][idx])["leverage"],
                             leverage_core=E.radius_covariance(Zc["fat"][idx], Zc["skin"][idx])["leverage"]))
    comp = pd.DataFrame(rows)
    Zf_, Zc_ = c.Z(wf), c.Z(core)
    rf = E.radius_covariance(Zf_["fat"][idx], Zf_["skin"][idx])
    rc = E.radius_covariance(Zc_["fat"][idx], Zc_["skin"][idx])
    rng = np.random.default_rng(M.SEED)
    nperm = 20 if chk else 2000
    nulls = {}
    for tag, Z_ in (("14gene", Zf_), ("core", Zc_)):
        rng = np.random.default_rng(M.SEED)  # same permutations for both waveforms
        zf, zs = Z_["fat"][idx], Z_["skin"][idx]
        nulls[tag] = np.array([E.radius_covariance(zf, zs[rng.permutation(len(zs))])["leverage"] for _ in range(nperm)])
    if chk:
        assert len(comp) == 6 and comp.leverage_core.notna().all()
        print("check ok")
        return
    orig = pd.read_csv(OUT.parents[1] / "P1_person_loudness/output/tests.csv").iloc[0]
    ocomp = pd.read_csv(OUT.parents[1] / "P1_person_loudness/output/companions.csv")
    assert np.allclose(comp.leverage_14gene, ocomp.leverage, atol=1e-9), "companions not reproduced"
    assert np.isclose(rf["leverage"], orig.leverage, atol=1e-9)
    pp = float((1 + np.sum(nulls["14gene"] >= rf["leverage"])) / (nperm + 1))
    if not np.isclose(pp, orig.p_perm, atol=1e-9):
        raise SystemExit(f"A331 permutation p not reproduced: stop before reading the core loop: {pp} vs {orig.p_perm}")
    print("14-gene reproduction OK", flush=True)
    out = []
    for tag, r, col in (("14gene", rf, "leverage_14gene"), ("core", rc, "leverage_core")):
        mt, rs = (comp[comp.kind == k][col] for k in ("matched", "responsive"))
        p = float((1 + np.sum(nulls[tag] >= r["leverage"])) / (nperm + 1))
        out.append(dict(waveform=tag, genes=len(wf.genes) if tag == "14gene" else len(core.genes), persons=r["n"], leverage=r["leverage"],
                        p_perm=p, matched_median=float(mt.median()), matched_q95=float(mt.quantile(0.95)),
                        pct_matched=float((mt < r["leverage"]).mean() * 100), responsive_median=float(rs.median()),
                        responsive_q95=float(rs.quantile(0.95)), pct_responsive=float((rs < r["leverage"]).mean() * 100)))
        out[-1]["holds"] = bool(r["leverage"] > out[-1]["matched_q95"] and p < 0.05)
    T = pd.DataFrame(out)
    T["original_leverage"] = [orig.leverage, np.nan]
    T["original_matched_q95"] = [orig.matched_q95, np.nan]
    T["original_responsive_q95"] = [orig.responsive_q95, np.nan]
    T.to_csv(OUT / "a335_1c_tests.csv", index=False)
    comp.to_csv(OUT / "a335_1c_companions.csv", index=False)
    pd.set_option("display.width", 200)
    print(T.round(4).T.to_string())


if __name__ == "__main__":
    main()
