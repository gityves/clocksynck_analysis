"""Adult power calibration of the planned A241 (fetal adrenal as host), before registration.

Read from the fetal data: the sample sheet only (`obs` of Cao 2020 `pseudobulk_organ.h5ad`: donor, compartment,
n_cells, depth). No fetal clock-gene count is read and no gate or coherence is run on fetal data.

Fetal design: pseudobulks with >= 1000 cells; the fetuses that have an Adrenal; their other organs. Each such fetus is
stood for by one distinct GTEx donor (death at 08-18 h, any Hardy class) having the adrenal and the stand-in organs of
that fetus (rarest sets first, `feasibility_fetal.draw`); other organs are masked. A fetus whose organ set is
carried by < 5 pool donors loses its rarest stand-in organ(s) until >= 5 do (once, in setup). Placenta, Eye, Thymus
and Spleen have no GTEx stand-in and are left out. Every stand-in pseudobulk is thinned binomially to the depth of the
fetal pseudobulk it stands for (`thinned`; `unthinned` otherwise).

Readout: `Dataset.gate(targets=other organs, hosts=["Adrenal"], min_shared=m)`, m in {5, 4, 3}; pass =
p_pairing < 0.05 and p_acrophase < 0.05; a ValueError = not estimable (NaN, not passed).
Controlled access: GTEx time of death is used in memory only; outputs are aggregates, without donor identifiers.

    python feasibility_a241.py [--n-sub 200] [--jobs 6] [--check]
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[2] / "scouts_2026-11-10/feasibility_fetal/scripts"))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402
from feasibility_fetal import GENES, ci, draw  # noqa: E402
from feasibility_gonad import thin  # noqa: E402

OUT = HERE.parents[1] / "output"
SHEET = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/natexp/cao2020/genome/pseudobulk_organ.h5ad"
STAND_IN = {"Adrenal": "Adrenal Gland", "Lung": "Lung", "Heart": "Heart - Left Ventricle", "Liver": "Liver",
            "Cerebellum": "Brain - Cerebellum", "Cerebrum": "Brain - Frontal Cortex (BA9)",
            "Muscle": "Muscle - Skeletal"}
ORGANS = ["Lung", "Placenta", "Heart", "Liver", "Cerebellum", "Cerebrum", "Muscle", "Eye", "Spleen", "Thymus"]
N_PERM, MIN_FULL, WF = 500, 5, cs.gtex()
POOL, FETAL, SETS = {}, None, []  # donor -> GTEx table; fetal rows (donor, organ, depth); fetal sets after masking


def setup():
    """Fetal sheet, GTEx pool, masking. Returns the design table."""
    global POOL, FETAL, SETS
    f = ad.read_h5ad(SHEET, backed="r").obs.astype({"donor": str, "compartment": str})
    f = f[f.n_cells >= 1000]
    f = f[f.donor.isin(f[f.compartment == "Adrenal"].donor)]
    pb, d = gtex.load(list(STAND_IN.values()), hardy=(0, 1, 2, 3, 4), depth="library_size")
    pb = pb[pb.donor.isin(set(d.index[in_window(d.tod.to_numpy())]))]
    pb = pb[pb.donor.isin(pb[pb.compartment == "Adrenal Gland"].donor)]
    POOL = {k: g.reset_index(drop=True) for k, g in pb.groupby("donor")}
    have = {k: set(g.compartment) for k, g in POOL.items()}
    n_pool = {o: sum(t in h for h in have.values()) for o, t in STAND_IN.items()}
    f = f.assign(tissue=f.compartment.map(STAND_IN))
    FETAL = f[f.tissue.notna()].reset_index(drop=True)
    masked = {o: 0 for o in STAND_IN}
    for k, g in FETAL.groupby("donor"):
        keep = set(g.tissue)
        while sum(keep <= h for h in have.values()) < MIN_FULL and len(keep) > 1:
            drop = min(keep - {"Adrenal Gland"}, key=lambda t: sum(t in h for h in have.values()))
            keep.discard(drop)
            masked[[o for o, t in STAND_IN.items() if t == drop][0]] += 1
            FETAL = FETAL.drop(FETAL.index[(FETAL.donor == k) & (FETAL.tissue == drop)])
    FETAL = FETAL.reset_index(drop=True)
    SETS = [set(g.tissue) for _, g in FETAL.groupby("donor", sort=True)]
    rows = []
    for o in ORGANS:
        s = f[f.compartment == o]
        rows.append(dict(organ=o, fetuses_with_adrenal=len(s), median_fetal_depth=s.depth.median() if len(s) else np.nan,
                         gtex_tissue=STAND_IN.get(o, "none"), pool_donors=n_pool.get(o, 0) if o in STAND_IN else 0,
                         fetuses_masked=masked.get(o, 0), fetuses_used=int((FETAL.compartment == o).sum())))
    return pd.DataFrame(rows)


def one(i):
    rng = np.random.default_rng([241, i])
    fet = [g.set_index("tissue") for _, g in FETAL.groupby("donor", sort=True)]
    pairs = draw(SETS, rng, POOL)
    tabs = {"thinned": [], "unthinned": []}
    for j, donor in pairs:
        t = POOL[donor]
        t = t[t.compartment.isin(SETS[j])].copy()
        t["donor"] = f"s{j}"
        tabs["unthinned"].append(t)
        tabs["thinned"].append(thin(t, t.compartment.map(fet[j].depth).to_numpy(float), rng))
    tg = [o for o in STAND_IN if o != "Adrenal"]
    out = []
    for name, ts in tabs.items():
        tab = pd.concat(ts, ignore_index=True)
        tab["compartment"] = tab.compartment.map({t: o for o, t in STAND_IN.items()})
        tab = tab[tab.groupby("compartment").donor.transform("nunique") >= 3]
        for m in (5, 4, 3):
            r = dict(variant=f"{name}_ms{m}", sub=i, n_pairs=np.nan, clock_share=np.nan, p_pairing=np.nan,
                     p_acrophase=np.nan)
            try:
                g = cs.Dataset(tab, anchors=None, waveform=WF, min_donors=3, name="a241").gate(
                    n_perm=N_PERM, seed=241 + i, targets=[o for o in tg if o in set(tab.compartment)],
                    hosts=["Adrenal"], min_shared=m)
                r.update(n_pairs=g.n_pairs, clock_share=g.clock_share, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase)
            except ValueError:
                pass
            out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0, sfx = time.time(), "_check" if a.check else ""
    limit_threads(1)
    OUT.mkdir(exist_ok=True)
    design = setup()
    design.to_csv(OUT / f"design{sfx}.csv", index=False)
    print(design.to_string(index=False), flush=True)
    n = 3 if a.check else a.n_sub
    res = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(one, range(n))):
            res.extend(r)
            if (k + 1) % 10 == 0:
                print(f"{k + 1}/{n} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res).sort_values(["variant", "sub"])
    res["passed"] = (res.p_pairing < 0.05) & (res.p_acrophase < 0.05)
    res.to_csv(OUT / f"subsamples{sfx}.csv", index=False)
    rows = []
    for v, r in res.groupby("variant", sort=False):
        k, m = int(r.passed.sum()), len(r)
        lo, hi = ci(k, m)
        rows.append(dict(variant=v, n_subsamples=m, pass_rate=k / m, ci_lo=lo, ci_hi=hi,
                         share_p_pairing=(r.p_pairing < 0.05).mean(), share_p_acrophase=(r.p_acrophase < 0.05).mean(),
                         estimable=r.p_pairing.notna().mean(), median_n_pairs=r.n_pairs.median()))
    rates = pd.DataFrame(rows)
    rates.to_csv(OUT / f"passrates{sfx}.csv", index=False)
    print(rates.round(3).to_string(index=False))
    print(f"runtime {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
