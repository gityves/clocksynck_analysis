"""Scout for a second misalignment design in blood (G2 (iii)): can the genome-wide rhythmic program of blood be placed on
central time or on the displaced sleep schedule, with the A261 machinery, in GSE48113 (Archer 2014, forced desynchrony)
and GSE122725 (Kervezee 2018 / 2019, simulated night shift, PBMC)?

Blood has no coherent core clock for this project (expected exception), so no host-referenced readout: the question is
only which time the rhythmic program follows. A scout reads designs, an external reference and noise; it does not
compute the program loadings of the test data on their real labels.

Reference (central time, no sleep, no meal rhythm): GSE39445 (Moller-Levet 2013), sleep-extension arm under extended
wakefulness, samples every 3 h with a melatonin-referenced circadian phase; per gene, cosinor on that phase with subject
intercepts and a linear term in hours awake (the wake drive). Reference-rhythmic: p < 0.01.
Program loading (time series analogue of A261): per subject and condition, the gene-centred values of the program's
genes are regressed on A_g (cos, sin)(w (t - phi_g)); z = a + i b, so a shift s and an excursion ratio r give
z = r e^{i w s}. Hypotheses (A261 rule, Hotelling of z_mis - r e^{i w s} z_al, level 0.01): central (s_c, 1),
behaviour (s_b, 1), lost (z_mis = 0).
Power: simulated programs (reference amplitudes x attenuation, phases with their standard errors) on the real design
(subjects, sample times), Gaussian noise per gene with the gene's SD within subject x condition (rhythm included, so
conservative); 200 programs per cell.

    python -I scout_blood.py
Outputs (output/): design.csv, reference_genes.csv, noise_summary.csv, power.csv.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
W = 2 * np.pi / 24
SEED, NREP, ALPHA = 1117, 200, 0.01
SIZES, ATTEN = [20, 50, 200, "all"], [1.0, 0.5, 0.25]


def series(path: Path):
    """Sample characteristics (DataFrame) and the expression table of a GEO series matrix."""
    meta, chars = {}, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                meta["gsm"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                break
        x = pd.read_csv(f, sep="\t", index_col=0, comment="!", low_memory=False)
    for c in chars:
        meta[c[0].split(":")[0].strip()] = [v.split(":", 1)[1].strip() if ":" in v else np.nan for v in c]
    x.index = x.index.astype(str)
    return pd.DataFrame(meta), x.apply(pd.to_numeric, errors="coerce")


def gpl15331():
    t = pd.read_csv(DATA / "blood_timed/GSE39445/GPL15331_table.txt", sep="\t", skiprows=21, dtype=str,
                    usecols=["ID", "GENE_SYMBOL"]).dropna(subset=["ID"])
    t = t[~t.ID.str.startswith("!")]
    sym = t.set_index("ID")["GENE_SYMBOL"].dropna()
    return sym[sym != ""]


def gpl24324():
    rx = re.compile(r"\(([^()\s]+)\),? (?:transcript variant [^,]+, )?(?:mRNA|non-coding RNA|long non-coding RNA)")
    out, on = {}, False
    with gzip.open(DATA / "kervezee2019_gse122725/GSE122725_family.soft.gz", "rt") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                on = True
                next(f)
                continue
            if line.startswith("!platform_table_end"):
                break
            if on:
                p = line.rstrip("\n").split("\t")
                if p[7] == "main":
                    m = rx.search(p[9].split(" /// ")[0])
                    if m:
                        out[p[0]] = m.group(1)
    return pd.Series(out)


def genes(x: pd.DataFrame, sym: pd.Series) -> pd.DataFrame:
    s = sym.reindex(x.index)
    return x[s.notna().to_numpy()].groupby(s.dropna().to_numpy()).mean()


def hours(s: str) -> float:
    h, m = s.split(":")
    return int(h) + int(m) / 60


def load_archer():
    smp, x = series(DATA / "archer2014_gse48113/GSE48113_series_matrix.txt.gz")
    smp = smp.assign(cond=np.where(smp["day sample was taken in protocol"] == "1", "al", "mis"),
                     t=[hours(v) if isinstance(v, str) and ":" in v else np.nan for v in smp["time sample taken"]])
    return smp.dropna(subset=["t"]), genes(x, gpl15331())


def load_kervezee():
    smp, x = series(DATA / "kervezee2019_gse122725/GSE122725_series_matrix.txt.gz")
    smp = smp[smp.experimentalgroup == "Control"]  # dim light; the bright-light arm (3 subjects) moves the clock
    smp = smp.assign(cond=np.where(smp.condition == "DayOrientedSchedule", "al", "mis"),
                     t=smp.relclocktime.astype(float))
    return smp, genes(x, gpl24324())


def reference(sym):
    smp, x = series(DATA / "blood_timed/GSE39445/GSE39445_series_matrix.txt.gz")
    smp = smp[smp.sleepprotocol == "Sleep Extension"].dropna(subset=["circadianphase"])
    y = genes(x[smp.gsm], sym).dropna()
    y = y[y.mean(1) > y.mean(1).quantile(0.2)]  # drop the lowest-intensity fifth
    ph, aw = smp.circadianphase.astype(float).to_numpy(), smp.hoursawake.astype(float).to_numpy()
    Dm = pd.get_dummies(smp.subject).to_numpy(float)
    D = np.c_[Dm, aw, np.cos(W * ph), np.sin(W * ph)]
    Y = y.to_numpy().T
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    res = Y - D @ B
    dfr = len(Y) - D.shape[1]
    s2 = (res ** 2).sum(0) / dfr
    D0 = D[:, :-2]
    rss0 = ((Y - D0 @ np.linalg.lstsq(D0, Y, rcond=None)[0]) ** 2).sum(0)
    p = stats.f.sf(((rss0 - (res ** 2).sum(0)) / 2) / s2, 2, dfr)
    XtXi = np.linalg.pinv(D.T @ D)
    a, b = B[-2], B[-1]
    return pd.DataFrame(dict(gene=y.index, amp=np.hypot(a, b), acro_h=np.mod(np.arctan2(b, a), 2 * np.pi) / W,
                             se=np.sqrt(s2 * (XtXi[-2, -2] + XtXi[-1, -1]) / 2), p=p,
                             subjects=smp.subject.nunique(), samples=len(smp)))


def loading(Y, t, amp, phi):
    """Y samples x genes (centred within the unit), t sample times: z = a + i b of Y ~ a X1 + b X2."""
    th = W * (t[:, None] - phi[None, :])
    X = np.c_[(amp * np.cos(th)).ravel(), (amp * np.sin(th)).ravel()]
    c, *_ = np.linalg.lstsq(X, Y.ravel(), rcond=None)
    return c[0] + 1j * c[1]


def classify(zA, zM, hyp):
    n, keep = len(zA), []
    for h, (s, r) in hyp.items():
        dd = zM - r * np.exp(1j * W * s) * zA
        X = np.c_[dd.real, dd.imag]
        m, S = X.mean(0), np.cov(X.T)
        T2 = n * m @ np.linalg.solve(S, m)
        if stats.f.sf((n - 2) / (2 * (n - 1)) * T2, 2, n - 2) >= ALPHA:
            keep.append(h)
    return keep[0] if len(keep) == 1 else ("mixed" if not keep else "unresolved")


def power(name, smp, sd, ref, hyp, rng):
    """Simulated programs on the design of `smp` (subjects with both conditions)."""
    units = {}
    for (s, c), g in smp.groupby(["subject", "cond"]):
        units.setdefault(s, {})[c] = g.t.to_numpy(float)
    units = {s: u for s, u in units.items() if len(u) == 2 and min(len(v) for v in u.values()) >= 3}
    rhy = ref[(ref.p < 0.01) & ref.gene.isin(sd.index)]
    rows = []
    for G in SIZES:
        for att in ATTEN:
            for truth, (s_t, r_t) in hyp.items():
                calls = []
                for _ in range(NREP):
                    g = rhy if G == "all" else rhy.iloc[rng.choice(len(rhy), G, replace=False)]
                    amp, sdg = g.amp.to_numpy(), sd.loc[g.gene].to_numpy()
                    phi = g.acro_h.to_numpy() + rng.normal(0, 1, len(g)) * g.se.to_numpy() / amp / W
                    a_t = np.maximum(amp + rng.normal(0, 1, len(g)) * g.se.to_numpy(), 0) * att
                    zA, zM = [], []
                    for u in units.values():
                        for c, out in (("al", zA), ("mis", zM)):
                            t = u[c]
                            sh, rr = (0.0, 1.0) if c == "al" else (s_t, r_t)
                            Y = rr * a_t * np.cos(W * (t[:, None] - phi - sh)) + rng.normal(0, 1, (len(t), len(g))) * sdg
                            out.append(loading(Y - Y.mean(0), t, amp, g.acro_h.to_numpy()))
                    calls.append(classify(np.array(zA), np.array(zM), hyp))
                c = pd.Series(calls).value_counts(normalize=True)
                rows.append(dict(design=name, G=len(rhy) if G == "all" else G, attenuation=att, truth=truth,
                                 units=len(units), correct=float(c.get(truth, 0.0)),
                                 **{f"call_{k}": float(v) for k, v in c.items()}))
                print(rows[-1], flush=True)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    design, noise, prows = [], [], []
    ref = reference(gpl15331())
    ref.to_csv(OUT / "reference_genes.csv", index=False)
    print("reference", len(ref), "genes,", int((ref.p < 0.01).sum()), "at p < 0.01", flush=True)
    for name, loader, hyp in (
            ("GSE48113 forced desynchrony", load_archer, {"central": (0.0, 1.0), "behaviour": (12.0, 1.0)}),
            ("GSE122725 night shift (dim-light control, = GSE107537)", load_kervezee, {"central": (0.0, 1.0), "behaviour": (10.0, 1.0)})):
        smp, y = loader()
        smp = smp[smp.gsm.isin(y.columns)]
        for c, g in smp.groupby("cond"):
            design.append(dict(design=name, cond=c, subjects=g.subject.nunique(), samples=len(g),
                               per_subject_median=float(g.groupby("subject").size().median()),
                               t_min=float(g.t.min()), t_max=float(g.t.max())))
        y = y[smp.gsm]
        yc = y.T.groupby([smp.subject.to_numpy(), smp.cond.to_numpy()]).transform(lambda v: v - v.mean())
        sd = yc.std().rename("sd")  # within subject x condition, rhythm included
        com = ref[ref.gene.isin(sd.index)]
        noise.append(dict(design=name, genes=len(y), reference_on_array=len(com),
                          reference_rhythmic_on_array=int((com.p < 0.01).sum()),
                          median_sd=float(sd.median()),
                          median_ref_amp_rhythmic=float(com[com.p < 0.01].amp.median()),
                          median_sd_rhythmic=float(sd.reindex(com[com.p < 0.01].gene).median())))
        print(noise[-1], flush=True)
        prows += power(name, smp, sd, ref, hyp | {"lost": (0.0, 0.0)}, rng)
    pd.DataFrame(design).to_csv(OUT / "design.csv", index=False)
    pd.DataFrame(noise).to_csv(OUT / "noise_summary.csv", index=False)
    pd.DataFrame(prows).to_csv(OUT / "power.csv", index=False)


if __name__ == "__main__":
    main()
