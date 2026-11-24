"""A365 (registered): do living persons differ in a clock amplitude shared by epidermis and dermis, measured from a
timed day? del Olmo et al. 2022, GSE205155 (A37 download): 11 subjects x 2 layers x 7 times (08:00 to 08:00 next day,
every 4 h), Agilent arrays, the submitters' normalised log2 values (series matrix), gene symbols from GPL14550.

One projector per gene set over all 154 arrays (GTEx not-sun-exposed skin waveform, log2 x ln 2, as A255 / A364);
states centred within subject x layer; per subject and layer the first-harmonic amplitude |beta| by least squares of
z_t = c + beta exp(i 2 pi t / 24) over the 7 times (independent of the person's phase); r = Pearson correlation over
subjects of log |beta| in epidermis and dermis. Core loop primary, 14 genes reported. Null: subjects permuted (10,000).
Companions: 200 expression-matched and 200 responsive pseudo-clock waveforms (seeds [365, 1] / [365, 2]). Core
readout: the 364 three-gene drops (CLAUDE.md gene-drop rule).

    python a365_delolmo_volume.py --stage power   # noise per layer, simulation, companions; no clock r
    python a365_delolmo_volume.py --stage score   # registered test (after the power note)
Outputs (`../output/`): a365_power.csv, a365_companions.csv, a365_sim.csv | a365_tests.csv, a365_drops.csv.
"""
import argparse
import gzip
import itertools
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import projector, pseudo_waveforms, states  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
from common.paths import DATA  # noqa: E402

SRC = DATA / "skin_paired" / "GSE205155"
OUT = HERE.parents[1] / "output"
WF = cs.gtex("Skin - Not Sun Exposed (Suprapubic)")
TRIO = ["PER1", "PER2", "NR1D1"]
CORE = WF.drop(TRIO)
CORE_IDX = [k for k, g in enumerate(WF.genes) if g in CORE.genes]
N_SETS, N_PERM, N_SIM = 200, 10_000, 2000
SIGMA_A = 0.352   # GTEx person spread of log amplitude (A331 power stage)
LAYERS = ("epidermis", "dermis")


def load():
    """(design: subject, layer, t in hours from the first 08:00; expression: arrays x gene symbols, natural log)."""
    meta, rows, head, in_tab = {}, [], None, False
    with gzip.open(SRC / "GSE205155_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_tab = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            p = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            if in_tab:
                if head is None:
                    head = p
                else:
                    rows.append(p)
            elif p[0] == "!Sample_title":
                meta["title"] = p[1:]
    x = pd.DataFrame([r[1:] for r in rows], index=[r[0] for r in rows], columns=head[1:]).astype(float)
    with open(SRC / "GPL14550_table.txt") as f:
        start = next(k for k, line in enumerate(f) if line.startswith("ID\t"))
    gpl = pd.read_csv(SRC / "GPL14550_table.txt", sep="\t", skiprows=start, usecols=["ID", "GENE_SYMBOL"],
                      comment="!", dtype=str)
    sym = gpl.set_index("ID").GENE_SYMBOL.dropna().astype(str).map(canonical)
    x = x[x.index.isin(sym.index)]
    x = x.groupby(sym.reindex(x.index).to_numpy()).mean().T * np.log(2.0)
    t = pd.Series(meta["title"]).str.extract(r"(\w+)_(\d+)h_subject(\d+)")
    d = pd.DataFrame(dict(layer=t[0].to_numpy(), t=t[1].astype(float).to_numpy() - 8.0, subject=t[2].to_numpy()))
    assert len(d) == 154 and set(d.layer) == set(LAYERS)
    x.index = range(len(d))
    return d, x


def site_states(x, d, w):
    g = list(w.genes)
    L = projector(x[g], w, counts=False)
    b, _ = states(x[g], L, w)
    z = b[:, 0] + 1j * b[:, 1]
    key = (d.subject + d.layer).to_numpy()
    return z - pd.Series(z).groupby(key).transform("mean").to_numpy()


def design(t):
    return np.column_stack([np.ones(len(t)), np.exp(1j * 2 * np.pi * np.asarray(t, float) / 24)])


def amplitudes(z, d):
    """|beta| per subject and layer (least squares on [1, exp(i w t)]); per-sample noise variance per layer."""
    out, res = {}, {lay: [] for lay in LAYERS}
    for (s, lay), i in d.groupby(["subject", "layer"]).groups.items():
        i = np.asarray(i)
        X = design(d.t.to_numpy()[i])
        coef = np.linalg.lstsq(X, z[i], rcond=None)[0]
        out[(s, lay)] = abs(coef[1])
        r = z[i] - X @ coef
        res[lay].append(float(np.sum(np.abs(r) ** 2) / (len(i) - 2)))
    a = pd.Series(out).unstack()[list(LAYERS)]
    return a, {lay: float(np.mean(v)) for lay, v in res.items()}


def r_of(a):
    la = np.log(a.to_numpy(float))
    return float(np.corrcoef(la[:, 0], la[:, 1])[0, 1])


def core_of(w, j):
    return Waveform(tuple(w.genes[k] for k in CORE_IDX), CORE.amp, CORE.psi, None, None, f"core{j}")


def simulate(m, s2, n, sigma, t, rng):
    X = design(t)
    A_ = np.exp(rng.normal(0, sigma, n))
    ph = rng.uniform(0, 2 * np.pi, n)
    la = []
    for lay in LAYERS:
        z = m[lay] * A_[:, None] * np.exp(1j * (2 * np.pi * t[None, :] / 24 + ph[:, None]))
        z = z + np.sqrt(s2[lay] / 2) * (rng.normal(size=z.shape) + 1j * rng.normal(size=z.shape))
        coef = np.linalg.lstsq(X, z.T, rcond=None)[0]
        la.append(np.log(np.abs(coef[1])))
    return float(np.corrcoef(la[0], la[1])[0, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    d, x = load()
    n = d.subject.nunique()
    t = np.sort(d.t.unique())
    missing = [g for g in WF.genes if g not in x.columns]
    if missing:
        raise SystemExit(f"clock genes missing: {missing}")
    lm = x.mean(0)
    lm = lm[~lm.index.duplicated()]
    print(f"{len(d)} arrays, {n} subjects, {x.shape[1]} genes, times {t.tolist()} ({time.time() - t0:.0f}s)", flush=True)

    if a.stage == "power":
        rows = []
        for kind, tag in (("matched", 1), ("responsive", 2)):
            for j, w in enumerate(pseudo_waveforms(lm, N_SETS, np.random.default_rng([365, tag]), WF, kind=kind)):
                for gs, ww in (("14 genes", w), ("core", core_of(w, j))):
                    am, _ = amplitudes(site_states(x, d, ww), d)
                    rows.append(dict(kind=kind, set=j, gene_set=gs, r=r_of(am)))
        comp = pd.DataFrame(rows)
        comp.to_csv(OUT / "a365_companions.csv", index=False)
        rng = np.random.default_rng([365, 3])
        XtX = np.linalg.inv(design(t).conj().T @ design(t))
        v22 = float(np.real(XtX[1, 1]))
        res, sims = {}, {}
        for gs, w in (("core", CORE), ("14 genes", WF)):
            am, s2 = amplitudes(site_states(x, d, w), d)
            m = {lay: float(np.sqrt(max(np.mean(am[lay].to_numpy() ** 2) - s2[lay] * v22, 1e-12))) for lay in LAYERS}
            null = np.array([simulate(m, s2, n, 0.0, t, rng) for _ in range(N_SIM)])
            for lab, sg in (("gtex", SIGMA_A), ("half", SIGMA_A / 2)):
                sims[f"{gs}|{lab}"] = np.array([simulate(m, s2, n, sg, t, rng) for _ in range(N_SIM)])
            c = comp[comp.gene_set == gs]
            mq95 = float(c[c.kind == "matched"].r.quantile(0.95))
            thr = max(mq95, float(np.quantile(null, 0.95)))
            res.update({f"{gs}|noise_epidermis": s2["epidermis"], f"{gs}|noise_dermis": s2["dermis"],
                        f"{gs}|m_epidermis": m["epidermis"], f"{gs}|m_dermis": m["dermis"],
                        f"{gs}|null_sim_q95": float(np.quantile(null, 0.95)),
                        f"{gs}|matched_median": float(c[c.kind == "matched"].r.median()), f"{gs}|matched_q95": mq95,
                        f"{gs}|responsive_median": float(c[c.kind == "responsive"].r.median()),
                        f"{gs}|responsive_q95": float(c[c.kind == "responsive"].r.quantile(0.95)),
                        f"{gs}|expected_r_gtex": float(np.median(sims[f'{gs}|gtex'])),
                        f"{gs}|power_gtex": float(np.mean(sims[f"{gs}|gtex"] > thr)),
                        f"{gs}|power_half": float(np.mean(sims[f"{gs}|half"] > thr))})
        res["gate"] = bool(res["core|power_gtex"] >= 0.8)
        pd.DataFrame(sims).to_csv(OUT / "a365_sim.csv", index=False)
        pd.Series(res).to_csv(OUT / "a365_power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "a365_power.csv", index_col=0)["value"]
    if str(pw["gate"]) != "True":
        raise SystemExit("A365 not powered: not decidable at this design; not scored")
    comp = pd.read_csv(OUT / "a365_companions.csv")
    rng = np.random.default_rng([365, 4])
    rows = []
    for gs, w in (("core", CORE), ("14 genes", WF)):
        am, _ = amplitudes(site_states(x, d, w), d)
        r = r_of(am)
        la = np.log(am.to_numpy(float))
        null = np.array([np.corrcoef(la[:, 0], la[rng.permutation(n), 1])[0, 1] for _ in range(N_PERM)])
        c = comp[comp.gene_set == gs]
        mt, rs = c[c.kind == "matched"].r, c[c.kind == "responsive"].r
        row = dict(test="H365.1" if gs == "core" else "H365.1 (14 genes, reported)", subjects=n, r=r,
                   p_perm=float((1 + np.sum(null >= r)) / (1 + N_PERM)), matched_median=float(mt.median()),
                   matched_q95=float(mt.quantile(.95)), pct_matched=float((mt < r).mean() * 100),
                   responsive_median=float(rs.median()), responsive_q95=float(rs.quantile(.95)),
                   pct_responsive=float((rs < r).mean() * 100), expected_r_gtex=float(pw[f"{gs}|expected_r_gtex"]))
        row["passes"] = bool(r > row["matched_q95"] and row["p_perm"] < 0.05)
        rows.append(row)
    dr = []
    for c3 in itertools.combinations(WF.genes, 3):
        am, _ = amplitudes(site_states(x, d, WF.drop(list(c3))), d)
        dr.append(dict(drop="+".join(c3), trio=list(c3) == TRIO, r=r_of(am)))
    dr = pd.DataFrame(dr)
    dr.to_csv(OUT / "a365_drops.csv", index=False)
    tr = float(dr[dr.trio].r.iloc[0])
    rows[0].update(drops_median=float(dr.r.median()), drops_q10=float(dr.r.quantile(.10)),
                   trio_drop_pct=float((dr.r < tr).mean() * 100))
    rows[0]["reading"] = ("living persons share a clock amplitude across epidermis and dermis, beyond matched genes"
                          if rows[0]["passes"] else "not shown beyond matched genes in timed living skin")
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "a365_tests.csv", index=False)
    print(out.round(4).T.to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
