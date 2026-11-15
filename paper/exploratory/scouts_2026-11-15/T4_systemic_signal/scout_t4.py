"""Scout for T4 (name the systemic signal): can per-sample readers of glucocorticoid and temperature signalling be built
from genes fixed outside the human data, and can a fixed phase relation between a reader and the host state be told from
none in GTEx? No hypothesis is scored: no observed reader-host relation is computed on the real pairing.

Part 1, mouse (Manella 2021, GSE159135; timed cosinors of A260, `A260_mouse_calibration/output/truth_genes.csv.gz`).
  Candidates (prompt): glucocorticoid FKBP5, TSC22D3, KLF9; cold-inducible CIRBP, RBM3; heat shock HSPA1A, HSPA1B, HSPA8,
  HSPH1, HSP90AA1, DNAJB1. Rule of the prompt: a candidate enters a reader if it is rhythmic (q < 0.05) in liver without
  a local clock (hepatocyte Bmal1 knockout, ad libitum). Expected phases in control liver: glucocorticoid targets after
  the corticosterone peak at the light-dark transition (about ZT 12-18); cold-inducible genes in the cool rest phase
  (about ZT 0-10); heat-shock genes in the warm active phase (about ZT 12-24).
Part 2, GTEx power (fast deaths; the 15 host tissues of the cue line, `gtex_cues/_shared.TISSUES`).
  Host per donor: `host_states` from the clock genes of all 15 tissues (and without PER1). Reader per sample: mean over
  the reader's genes of the within-tissue z of log CPM. Statistic: per tissue, the complex regression of the centred
  reader on the centred host state, beta_k; consistency K = |sum_k beta_k| / sum_k |beta_k| over tissues (1: one phase
  relation in every tissue); p against 200 permutations of the reader across donors within tissue.
  Simulation: the real reader permuted across donors within tissue (real relation removed, marginals kept) plus a
  relation of fixed phase (random per draw, the same in every tissue) scaled to an R^2 on the host of 0.014 (the A27 GC
  score in fast deaths), 0.03, 0.05 and 0 (false-positive rate); 100 draws each. Brain looser: descriptive.
Controlled access: GTEx time of death is not used; outputs are aggregates.

    python scout_t4.py [--draws 100] [--jobs 6]
Outputs (output/): mouse_candidates.csv, power.csv.
"""
import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_3/scripts/gtex_cues"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE / "output"
SEED, N_PERM = 1118, 200
READERS = {"glucocorticoid": ("FKBP5", "TSC22D3", "KLF9"), "cold": ("CIRBP", "RBM3"),
           "heat": ("HSPA1A", "HSPA1B", "HSPA8", "HSPH1", "HSP90AA1", "DNAJB1")}
EXPECTED = {"glucocorticoid": (12, 18), "cold": (0, 10), "heat": (12, 24)}
R2 = (0.0, 0.014, 0.03, 0.05)
G = {}


def mouse():
    t = pd.read_csv(PAPER / "exploratory/A260_mouse_calibration/output/truth_genes.csv.gz")
    rows = []
    for r, genes in READERS.items():
        for g in genes:
            row = dict(reader=r, gene=g)
            for arm in ("Liver|AlbCre|AL", "Liver|BLKO|AL", "Liver|AlbCre|DF", "Liver|BLKO|DF"):
                x = t[(t.arm == arm) & (t.gene == g)]
                tag = arm.split("|", 1)[1].replace("|", "_")
                row |= {f"q_{tag}": float(x.q.iloc[0]) if len(x) else np.nan,
                        f"phase_{tag}": float(x.phase_h.iloc[0]) if len(x) else np.nan}
            lo, hi = EXPECTED[r]
            row["rule_rhythmic_without_clock"] = bool(row["q_BLKO_AL"] < 0.05)
            row["control_rhythmic_at_expected_phase"] = bool(row["q_AlbCre_AL"] < 0.05 and lo <= row["phase_AlbCre_AL"] % 24 < hi)
            rows.append(row)
    return pd.DataFrame(rows)


def consistency(reader: pd.DataFrame, host: pd.DataFrame, rng, n_perm):
    """K and its permutation p; reader: rows (donor, compartment, value)."""
    betas, data = [], []
    for k, g in reader.groupby("compartment"):
        g = g[g.donor.isin(host.index)]
        z = (host.loc[g.donor, "x"].to_numpy() + 1j * host.loc[g.donor, "y"].to_numpy())
        z = z - z.mean()
        y = g.value.to_numpy() - g.value.mean()
        data.append((y, z))
    def K(perm):
        b = np.array([np.sum((y[p] if p is not None else y) * np.conj(z)) / np.sum(np.abs(z) ** 2)
                      for (y, z), p in zip(data, perm)])
        return np.abs(b.sum()) / np.abs(b).sum()
    obs = K([None] * len(data))
    null = np.array([K([rng.permutation(len(y)) for y, _ in data]) for _ in range(n_perm)])
    return obs, (1 + np.sum(null >= obs)) / (1 + n_perm)


def draw(args):
    reader, r2, i = args
    rng = np.random.default_rng([SEED, i, int(r2 * 1000), list(READERS).index(reader)])
    host, R = G["host"], G["reader"][reader]
    rows = []
    phase = rng.uniform(0, 2 * np.pi)
    sim = []
    for k, g in R.groupby("compartment"):
        g = g[g.donor.isin(host.index)].copy()
        v = rng.permutation(g.value.to_numpy())
        z = host.loc[g.donor, "x"].to_numpy() + 1j * host.loc[g.donor, "y"].to_numpy()
        s = np.real(np.exp(-1j * phase) * (z - z.mean()))
        if r2 > 0:  # signal scaled so that its share of the reader variance is r2
            s = s / s.std() * v.std() * np.sqrt(r2 / (1 - r2))
            v = v + s
        g["value"] = v
        sim.append(g)
    K, p = consistency(pd.concat(sim), host, rng, N_PERM)
    return dict(reader=reader, r2=r2, draw=i, K=K, p=p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    m = mouse()
    m.to_csv(OUT / "mouse_candidates.csv", index=False)
    print(m.round(3).to_string(index=False), flush=True)
    genes = sorted({g for v in READERS.values() for g in v})
    pb, _ = gtex.load(tissues=list(S.TISSUES), hardy=(1, 2), genes=[*cs.gtex().genes, *genes], depth="library_size",
                      require_tod=False)
    G["host"] = S.host_states(pb[["donor", "compartment", "depth", *cs.gtex().genes]], S.TISSUES, name="fast")
    lg = np.log(pb[genes].to_numpy(float) / pb.depth.to_numpy(float)[:, None] * 1e6 + 0.5)
    z = pd.DataFrame(lg, columns=genes).groupby(pb.compartment.to_numpy()).transform(lambda x: (x - x.mean()) / x.std())
    G["reader"] = {r: pd.DataFrame(dict(donor=pb.donor, compartment=pb.compartment, value=z[list(g)].mean(1)))
                   for r, g in READERS.items()}
    print(f"host donors {len(G['host'])}; samples {len(pb)}", flush=True)
    tasks = [(r, r2, i) for r in READERS for r2 in R2 for i in range(a.draws)]
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        res = pd.DataFrame(p.map(draw, tasks))
    pw = res.groupby(["reader", "r2"]).agg(draws=("p", "size"), pass_rate=("p", lambda p: float((p < 0.05).mean())),
                                           median_K=("K", "median")).reset_index()
    pw.to_csv(OUT / "power.csv", index=False)
    print(pw.round(3).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
