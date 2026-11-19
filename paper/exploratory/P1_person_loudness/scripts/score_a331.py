"""A331 (registered, check P1): do living persons differ in a clock loudness shared by fat and skin (MuTHER)?

    python score_a331.py --stage power   # GTEx person spread, plane-level simulation, companion leverages; no clock radius
    python score_a331.py --stage score   # clock leverage, percentile, twin-pair bootstrap, zygosity (if the gate passed)
Outputs (`../output/`): power.csv, sim.csv, companions.csv | tests.csv. Aggregates only.
"""
import argparse
import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import elasticity as E  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, N_SIM, N_BOOT = 331, 1000, 2000
F_FAT, F_SKIN = 0.696, 0.408  # clock fractions of state variance (A255 gate shares and pairing 0.533)


def a255():
    spec = importlib.util.spec_from_file_location("a255", PAPER / "exploratory/A255_living_tissues/scripts/compute_a255.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def gtex_sigma():
    import compute_oscillator_gtex as A104
    from common import gtex
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS)
    pb = pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)
    comps = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb, anchors=comps, name="gtex_periphery")
    cv = []
    for hA, hB in E.split_hosts(ds, comps, 20, seed=[SEED, 1]):
        j = hA.index.intersection(hB.index)
        r = E.radius_covariance((hA.x + 1j * hA.y)[j].to_numpy(), (hB.x + 1j * hB.y)[j].to_numpy())
        cv.append(r["cv2"])
    return float(np.sqrt(max(np.median(cv), 0))), cv


def simulate(n, sigma, phases, rng):
    A = np.exp(rng.normal(0, sigma, n) - sigma ** 2 / 2)
    T = A * np.exp(1j * rng.choice(phases, n))
    T = T - T.mean()
    v = np.mean(np.abs(T) ** 2)
    out = []
    for f in (F_FAT, F_SKIN):
        s = np.sqrt(v * (1 - f) / f / 2)
        out.append(T + s * (rng.normal(size=n) + 1j * rng.normal(size=n)))
    return E.radius_covariance(*out)["leverage"]


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    m = a255()
    c, persons = m.load_muther(False)
    both = persons[(persons.has_fat == "True") & (persons.has_skin == "True")].person
    idx = np.array([i in set(both) for i in c.ids])
    if a.stage == "power":
        sigma, cv = gtex_sigma()
        Z = c.Z()
        zf, zs = Z["fat"][idx], Z["skin"][idx]
        rot = np.angle(np.nansum(zf * np.conj(zs)))
        phases = np.angle(zf / np.abs(zf).mean() + zs * np.exp(1j * rot) / np.abs(zs).mean())
        rng = np.random.default_rng(SEED)
        n = int(idx.sum())
        sims = {s: np.array([simulate(n, s, phases, rng) for _ in range(N_SIM)]) for s in (0.0, sigma / 2, sigma)}
        rows = []
        for kind, tag in (("matched", 3), ("responsive", 4)):
            for j, w in enumerate(c.pseudo(kind, m.N_SETS, tag)):
                Zp = c.Z(w)
                rows.append(dict(kind=kind, set=j, leverage=E.radius_covariance(Zp["fat"][idx], Zp["skin"][idx])["leverage"]))
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "companions.csv", index=False)
        mq95 = float(comp[comp.kind == "matched"].leverage.quantile(0.95))
        rq95 = float(comp[comp.kind == "responsive"].leverage.quantile(0.95))
        null_q95 = float(np.quantile(sims[0.0], 0.95))
        pw = dict(persons=n, gtex_cv2_median=float(np.median(cv)), sigma_A=sigma, null_sim_q95=null_q95,
                  matched_median=float(comp[comp.kind == "matched"].leverage.median()), matched_q95=mq95,
                  responsive_median=float(comp[comp.kind == "responsive"].leverage.median()), responsive_q95=rq95,
                  expected_leverage_sigma=float(np.median(sims[sigma])),
                  expected_leverage_half=float(np.median(sims[sigma / 2])),
                  sim_q05_sigma=float(np.quantile(sims[sigma], 0.05)),
                  power_vs_matched_q95=float(np.mean(sims[sigma] > mq95)),
                  power_vs_sim_null=float(np.mean(sims[sigma] > null_q95)),
                  power_half_vs_matched_q95=float(np.mean(sims[sigma / 2] > mq95)))
        pw["gate"] = bool(pw["power_vs_matched_q95"] >= 0.8)
        pd.Series(pw).to_csv(OUT / "power.csv", header=["value"])
        pd.DataFrame({f"sigma_{k:.3f}": v for k, v in sims.items()}).to_csv(OUT / "sim.csv", index=False)
        print(pd.Series(pw).to_string())
        return
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("power gate not passed: not decidable at this design; not scored")
    Z = c.Z()
    ids = np.array(c.ids)[idx]
    zf, zs = Z["fat"][idx], Z["skin"][idx]
    r = E.radius_covariance(zf, zs)
    rng = np.random.default_rng(SEED)
    null = np.array([E.radius_covariance(zf, zs[rng.permutation(len(zs))])["leverage"] for _ in range(2000)])
    comp = pd.read_csv(OUT / "companions.csv")
    mt, rs = comp[comp.kind == "matched"].leverage, comp[comp.kind == "responsive"].leverage
    # twin-pair bootstrap: a unit = a pair (both co-twins) or a singleton
    p = persons.set_index("person")
    unit = {}
    for i in ids:
        tw = p.cotwin.get(i)
        unit[i] = "|".join(sorted([i, tw])) if isinstance(tw, str) and tw in set(ids) else i
    u = pd.Series(unit)
    groups = {k: np.flatnonzero(np.isin(ids, v.index)) for k, v in u.groupby(u)}
    keys = list(groups)
    boot = []
    for _ in range(N_BOOT):
        pick = np.concatenate([groups[k] for k in rng.choice(keys, len(keys))])
        boot.append(E.radius_covariance(zf[pick], zs[pick])["leverage"])
    lo, hi = np.quantile(boot, [0.05, 0.95])
    t = dict(test="H331.1", persons=r["n"], leverage=r["leverage"], cv2=r["cv2"], lo90_twin_pairs=float(lo),
             hi90_twin_pairs=float(hi), p_perm=float((1 + np.sum(null >= r["leverage"])) / 2001),
             matched_median=float(mt.median()), matched_q95=float(mt.quantile(0.95)),
             pct_matched=float((mt < r["leverage"]).mean() * 100), responsive_median=float(rs.median()),
             responsive_q95=float(rs.quantile(0.95)), pct_responsive=float((rs < r["leverage"]).mean() * 100),
             sim_q05_sigma=float(pw["sim_q05_sigma"]), expected_leverage_sigma=float(pw["expected_leverage_sigma"]))
    t["verdict"] = ("PASS" if t["leverage"] > t["matched_q95"] and t["p_perm"] < 0.05 else
                    "specific to deaths" if t["leverage"] < t["sim_q05_sigma"] else
                    "not separable from a shared technical scale")
    rows = [t]
    zy = p.zygosity.reindex(ids)
    for z in ("MZ", "DZ"):
        k = (zy == z).to_numpy()
        if k.sum() > 30:
            rows.append(dict(test=f"descriptive {z}", persons=int(k.sum()), leverage=E.radius_covariance(zf[k], zs[k])["leverage"]))
    pd.DataFrame(rows).to_csv(OUT / "tests.csv", index=False)
    print(pd.DataFrame(rows).round(4).T.to_string())


if __name__ == "__main__":
    main()
