"""A345 (registered in `paper/docs/registered_experiments.md`): is the psoriatic lesional clock delay carried by
pre-mRNA? Tsoi 2019 (GSE121212), lesional vs non-lesional pairs; psoriasis primary, atopic dermatitis descriptive.

Stages:
- `power`: unspliced stand-in = each run's spliced clock counts thinned binomially to the run's genome-wide
  direct-intron / exon depth ratio (`qc.tsv`); 500 draws; per draw `clocksync.paired` delay p and gate, and the
  offset difference D (thinned vs spliced) with a patient bootstrap. No intron count of any gene is read.
- `check`: the scoring code with unspliced := spliced thinned to 10 % and the lesional / non-lesional labels swapped
  within a random half of the patients (no real intron count; expected offset near 0).
- `score`: the registered run. Spliced = `counts_exon`, unspliced = `counts_intron_direct` (`counts_intron` as
  sensitivity). Per readout: `clocksync.paired` (unit = patient, non-lesional -> lesional, GTEx not-sun-exposed skin
  waveform, 14 genes and core loop; gate, offset, one-sided delay p, 200 matched + 200 responsive companions on r);
  offset companions: the coherence-weighted delay |c| sin(arg c) on 200 matched and 200 responsive pseudo-clock
  waveforms (pre-data correction, registration note); D = offset_u - offset_s
  with 2,000 patient bootstraps (both readouts resampled jointly).
P345.0: spliced delay p < 0.05 and gate p < 0.05. H345.1: unspliced gate p < 0.05 and delay p < 0.05. H345.2: 90 %
interval of D inside +-1 h. H345.3: unspliced |c| sin(arg c) below the q05 of both companion families.

    python a345.py --stage power|check|score [--dataset tsoi2019_skin_trimmed]
Outputs (`../output/<stage>/`).
"""
import argparse
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
from clocksync.paired import ccorr, cohort_states, pseudo_waveforms  # noqa: E402

ROOT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/sccc_mapping")
OUT = HERE.parents[1] / "output"
SEED, N_PERM, N_BOOT, N_COMP, N_DRAW = 3450, 10000, 2000, 200, 500
HOURS = 24 / (2 * np.pi)
WF = cs.gtex("Skin - Not Sun Exposed (Suprapubic)")
WAVEFORMS = {"all": WF, "core": WF.drop(["PER1", "PER2", "NR1D1"])}
LEVELS = ("non-lesional", "lesional")
DISEASES = ("PSO", "AD")


def design(ds) -> pd.DataFrame:
    s = pd.read_csv(ROOT / ds / "samples.tsv", sep="\t")
    q = pd.read_csv(ROOT / ds / "qc.tsv", sep="\t").set_index("run")
    t = s.geo_title.str.extract(r"^(PSO|AD|CTRL)_(\d+)_(.*)$")
    d = pd.DataFrame(dict(run=s.run_accession, disease=t[0], patient=t[0] + "_" + t[1], condition=t[2]))
    d = d[d.disease.isin(DISEASES) & d.condition.isin(LEVELS)]
    both = d.groupby("patient").condition.nunique()
    d = d[d.patient.isin(both.index[both == 2])].reset_index(drop=True)
    d["ratio"] = (q.loc[d.run, "intron_direct_assigned_n"] / q.loc[d.run, "assigned_exon_n"]).to_numpy()
    return d


def load(ds, name) -> pd.DataFrame:
    x = pd.read_csv(ROOT / ds / f"counts_{name}.tsv.gz", sep="\t")
    x["symbol"] = x.symbol.astype(str).map(canonical)
    return x.drop_duplicates("symbol").set_index("symbol").drop(columns="gene_id")


def frame(d, X, genes) -> pd.DataFrame:
    """paired() input: one row per sample, unit, condition, depth and gene counts."""
    f = X.loc[genes, d.run].T.reset_index(drop=True)
    f.insert(0, "depth", X[d.run].sum().to_numpy(float))
    f.insert(0, "condition", d.condition.to_numpy())
    f.insert(0, "unit", d.patient.to_numpy())
    return f


def offset(f, wf, idx=None, signed=False) -> float:
    """Lesional minus non-lesional offset (h) of the pairs in f (rows of `idx` patients if given); `signed`: the
    coherence-weighted delay component |c| sin(arg c) instead (near 0 without a shared rhythm; negative = delay)."""
    case = f[f.condition == "lesional"].set_index("unit")
    ref = f[f.condition == "non-lesional"].set_index("unit").loc[case.index]
    if idx is not None:
        case, ref = case.loc[idx], ref.loc[idx]
    bc, vc, br, vr = cohort_states(case, ref, True, wf)
    c = ccorr(bc - bc.mean(0), vc, br - br.mean(0), vr)
    return float(np.imag(c)) if signed else float(np.angle(c) * HOURS)


def wrap(h):
    return (np.asarray(h) + 12) % 24 - 12


def interval(est, draws):
    x = wrap(np.asarray(draws) - est) + est
    return float(np.quantile(x, .05)), float(np.quantile(x, .95))


def readout(f, wf, seed, means=None, n_perm=N_PERM) -> dict:
    kw = dict(pseudo_means={"all": means}, pseudo_responsive=True) if means is not None else {}
    r = cs.paired(f, unit="unit", condition="condition", levels=LEVELS, waveform=wf, n_perm=n_perm, n_boot=1000,
                  gated_only=False, min_pairs=5, seed=seed, **kw)
    c = r.cohorts.iloc[0]
    out = dict(pairs=int(c.n_pairs), gate_p=float(c.gate_p), offset_h=float(c.offset_h), r=float(c.r),
               p_delay=float(r.pooled.get("p_offset_lower", np.nan)))
    return out


def comp_offsets(d, X, wf, means, rng) -> dict:
    out = {}
    for kind in ("matched", "responsive"):
        ws = pseudo_waveforms(means, N_COMP, rng, wf, kind)
        offs = []
        for w in ws:
            g = list(w.genes)
            offs.append(offset(frame(d, X, g), w, signed=True))
        out[f"{kind}_q05"], out[f"{kind}_q95"] = np.quantile(offs, [.05, .95])
    return out


def score(ds, check):
    out = OUT / ("check" if check else "score")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    d = design(ds)
    ex, itd, it = load(ds, "exon"), load(ds, "intron_direct"), load(ds, "intron")
    if check:
        th = pd.DataFrame(rng.binomial(ex.to_numpy(np.int64), 0.1), index=ex.index, columns=ex.columns)
        itd = it = th
        for p in rng.choice(d.patient.unique(), d.patient.nunique() // 2, replace=False):
            m = d.patient == p
            d.loc[m, "condition"] = d.loc[m, "condition"].map({"lesional": "non-lesional", "non-lesional": "lesional"})
    rows = []
    for dis in DISEASES:
        dd = d[d.disease == dis].reset_index(drop=True)
        means = np.log(ex[dd.run].div(ex[dd.run].sum(), axis=1).mean(1) + 1e-12)
        means = means[ex[dd.run].mean(1) >= 5]
        pts = dd.patient.unique()
        B = [rng.choice(pts, len(pts)) for _ in range(N_BOOT)]
        for w, wf in WAVEFORMS.items():
            g = list(wf.genes)
            fs, fu, fs2 = frame(dd, ex, g), frame(dd, itd, g), frame(dd, it, g)
            allg = list(dict.fromkeys(g + [x for x in means.index if x in itd.index]))  # companions need their genes
            row = dict(disease=dis, waveform=w)
            for k, X in (("s", ex), ("u", itd)):
                row.update({f"{k}_{a}": b for a, b in readout(frame(dd, X, allg), wf, SEED, means).items()})
            row.update({f"u_sub_{a}": b for a, b in readout(fs2, wf, SEED).items()})
            os_b = np.array([offset(fs, wf, b) for b in B])
            ou_b = np.array([offset(fu, wf, b) for b in B])
            row["D"] = float(wrap(row["u_offset_h"] - row["s_offset_h"]))
            row["D_lo"], row["D_hi"] = interval(row["D"], wrap(ou_b - os_b))
            row["u_lo"], row["u_hi"] = interval(row["u_offset_h"], ou_b)
            row["s_lo"], row["s_hi"] = interval(row["s_offset_h"], os_b)
            row["u_signed"] = offset(fu, wf, signed=True)
            row.update({f"u_comp_{k}": v for k, v in comp_offsets(dd, itd, wf, means, rng).items()})
            rows.append(row)
            print(dis, w, flush=True)
    t = pd.DataFrame(rows)
    t["P345_0"] = (t.s_p_delay < .05) & (t.s_gate_p < .05)
    t["H345_1"] = (t.u_p_delay < .05) & (t.u_gate_p < .05)
    t["H345_2"] = (t.D_lo > -1) & (t.D_hi < 1)
    t["H345_3"] = (t.u_signed < t.u_comp_matched_q05) & (t.u_signed < t.u_comp_responsive_q05)
    t.to_csv(out / "results.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.round(3).T.to_string())


def power(ds):
    out = OUT / "power"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng([SEED, 1])
    d = design(ds)
    ex = load(ds, "exon")
    rows = []
    dd = d[d.disease == "PSO"].reset_index(drop=True)
    pts = dd.patient.unique()
    for w, wf in WAVEFORMS.items():
        g = list(wf.genes)
        fs = frame(dd, ex, g)
        s_off = offset(fs, wf)
        for i in range(N_DRAW):
            th = rng.binomial(ex.loc[g, dd.run].to_numpy(np.int64), dd.ratio.to_numpy()[None, :])
            fu = fs.copy()
            fu[g] = th.T
            fu["depth"] = fs.depth.to_numpy() * dd.ratio.to_numpy()
            r = readout(fu, wf, SEED + i, n_perm=2000)
            B = [rng.choice(pts, len(pts)) for _ in range(200)]
            Db = wrap(np.array([offset(fu, wf, b) - offset(fs, wf, b) for b in B]))
            D = float(wrap(r["offset_h"] - s_off))
            lo, hi = interval(D, Db)
            rows.append(dict(waveform=w, draw=i, s_offset=s_off, u_offset=r["offset_h"], gate_p=r["gate_p"],
                             p_delay=r["p_delay"], D=D, D_halfwidth=(hi - lo) / 2))
            if (i + 1) % 50 == 0:
                print(w, i + 1, flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(out / "draws.csv", index=False)
    s = t.groupby("waveform").agg(s_offset=("s_offset", "first"), u_offset_median=("u_offset", "median"),
                                  u_offset_sd=("u_offset", "std"),
                                  h1_rate=("p_delay", lambda x: float(np.mean((x < .05) & (t.loc[x.index, "gate_p"] < .05)))),
                                  D_halfwidth_median=("D_halfwidth", "median")).reset_index()
    s["H345_1_decidable"], s["H345_2_decidable"] = s.h1_rate >= .8, s.D_halfwidth_median < 1
    s.to_csv(out / "summary.csv", index=False)
    print(s.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("power", "check", "score"), required=True)
    ap.add_argument("--dataset", default="tsoi2019_skin_trimmed")
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    if a.stage == "power":
        power(a.dataset)
    else:
        score(a.dataset, a.stage == "check")
    print(f"runtime {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
