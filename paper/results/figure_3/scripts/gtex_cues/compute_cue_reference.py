"""A81 and A81b: a GTEx cue reference, defined by ventilated vs fast deaths, and its calibration.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A81; A81b written after the A81 selection and evaluation). Estimand: per
gene and tissue, the timed cosinor loading in fast deaths (c_fast, p_fast) and the loading on the host clock state in
ventilated donors (c_int); lock = |c_int| / |c_fast| relative to the clock-plane genes of the tissue. Cue set = low lock,
locked set = high lock (clock outputs).
- A81: donors with a time of death split in halves (seed 81); selection halves define the sets (eligible: >= 4 rhythmic
  tissues, phase coherence R >= 0.7; lowest / highest lock tertile, top 40 by amplitude), evaluation halves score them:
  E1 (core-vs-set offsets per tissue, fast / ventilated) and H81.1 (the core-cue coupling drops in ventilated donors more
  than the core-locked coupling; Hardy labels permuted within tissue, 2,000). Expected: 3-gene sets, uninformative.
- A81b: rhythm from all fast deaths, lock from the ventilated selection half; 11 tissues (without whole blood);
  eligible nt >= 8, Fisher p < 1e-4, Rw >= 0.6; cue = below the median lock, locked = at or above, top 40 each; frame check
  (the 14 clock genes processed like the others must reproduce the fixed waveform; legacy: failed at 0.48, overridden post
  hoc), E1 on the evaluation halves. H81.1 is not confirmatory for A81b.
Hosts: `clocksync` `Dataset.fit(...).host`; core and set states: `clocksync.paired` (one projector, Poisson-corrected
complex correlation). Aggregates only (GTEx controlled access).

    python compute_cue_reference.py [--subsample N] [--jobs 2]
Outputs (paper/results/figure_3/output/gtex_cues/): a81_selection_genes.csv.gz, a81_sets.csv, a81_e1.csv,
a81_h1_verdict.csv, selection_genes_26.csv.gz, a81b_selection_genes.csv.gz, a81b_frame_check.csv, a81b_sets.csv,
a81b_e1.csv, **a81b_cue_reference.csv** (+ a81b_cue_reference.README.txt): the cue reference used by other figures.
Legacy: legacy/desync/validation/paired_damping/scripts/a81_gtex_cues.py, a81b_cues.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import importlib  # noqa: E402

P = importlib.import_module("clocksync.paired")  # the module (clocksync.paired is also the function)
from clocksync.genes import gene_columns  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _genes as G  # noqa: E402
from _genes import FOLDER, W  # noqa: E402

README = """a81b_cue_reference.csv: the GTEx cue reference of amendment A81b (GTEx v10, fast deaths for rhythm, ventilated donors
for lock; 11 tissues; aggregates over genes and tissues, no per-donor data).

columns
  gene      gene symbol (GTEx / HGNC)
  set       cue (lock below the median of the eligible genes: rhythm follows timing cues, not the local clock) or locked
            (lock at or above the median: clock-output-like)
  nt        tissues in which the gene is expressed (median CPM >= 5)
  fisher    Fisher combined p of the timed cosinor p-values over tissues
  Rw        weighted phase coherence |sum c| / sum |c| of the per-tissue complex loadings
  lock      median lock over tissues with p_fast < 0.05
  amp       amplitude |sum c| / nt (log CPM units per cosine)
  acro_h    acrophase (h, time of death = clock time, peak of cos(w (t - acro_h))), GTEx frame of the core waveform

frame check (a81b_frame_check.csv): {frame}

build a clocksync waveform (as the legacy A81.waveforms):
    ref = pd.read_csv("a81b_cue_reference.csv", index_col=0)
    wf = {{k: cs.waveform.Waveform(tuple(x.index), x.amp.to_numpy(), x.acro_h.to_numpy() * 2 * np.pi / 24, None, None, k)
          for k, x in ref.groupby("set")}}
(genes absent from a dataset are dropped by `Waveform.subset`). A constant frame shift (the frame check gives its size)
cancels in two-sample comparisons on the same waveform (paired designs); it enters within-sample offsets.
"""


def eval_tables(tissues, H, wfs) -> dict:
    """Per tissue: counts of core and set genes, depth, Hardy group and time of death of the evaluation donors."""
    ev = H["fast"]["eval"] | H["vent"]["eval"]
    out = {}
    for tis in tissues:
        X, g, ob = G.load_rows(tis, (0, 1, 2), ev)
        names = {}
        for wf in wfs.values():
            idx = gene_columns(g, wf.genes)
            for n, c in zip(wf.genes, idx):
                if c >= 0:
                    names[n] = c
        d = pd.DataFrame(np.round(X[:, list(names.values())] * ob.library_size.to_numpy(float)[:, None] / 1e6),
                         columns=list(names))
        d["depth"] = ob.library_size.to_numpy(float)
        d["vent"] = (ob.hardy == 0).to_numpy()
        d["tod"] = ob.tod.to_numpy(float)
        out[tis] = d
    return out


def set_states(d, wfs) -> dict:
    st = {}
    for k, wf in wfs.items():
        L = P.projector(d, wf, counts=True)
        st[k] = P.states(d, L, wf, depth=d.depth.to_numpy(float))
    return st


def corr(st, m, k):
    (bc, vc), (bs, vs) = st["core"], st[k]
    return P.ccorr(bc[m], vc[m], bs[m], vs[m])


def evaluate(tables, wfs, n_perm, seed=81):
    rows, states = [], {}
    for tis, d in tables.items():
        st = set_states(d, wfs)
        states[tis] = (st, d.vent.to_numpy())
        for grp, m in (("fast", ~d.vent.to_numpy()), ("vent", d.vent.to_numpy())):
            row = dict(tissue=tis, group=grp, n=int(m.sum()))
            for k in ("cue", "locked"):
                r = corr(st, m, k)
                row[f"{k}_offset_h"], row[f"{k}_corr"] = float(np.angle(r) * P.HOURS), float(abs(r))
            t = d.tod.to_numpy()[m]
            for k in ("core", "cue", "locked"):
                b = st[k][0][m]
                c = b[:, 0] + 1j * b[:, 1]
                X = np.column_stack([np.ones(m.sum()), np.exp(1j * W * t), np.exp(-1j * W * t)])
                row[f"{k}_timed_acro_h"] = float(np.angle(np.linalg.lstsq(X, c, rcond=None)[0][1]) / W % 24)
            rows.append(row)
    e1 = pd.DataFrame(rows)

    def stat(labels):
        drop = {k: np.mean([abs(corr(st, ~labels[t], k)) - abs(corr(st, labels[t], k)) for t, (st, _) in
                            states.items()]) for k in ("cue", "locked")}
        return drop["cue"] - drop["locked"], drop

    obs, drops = stat({t: v for t, (_, v) in states.items()})
    rng = np.random.default_rng(seed)
    null = np.array([stat({t: rng.permutation(v) for t, (_, v) in states.items()})[0] for _ in range(n_perm)])
    v = dict(D=obs, drop_cue=drops["cue"], drop_locked=drops["locked"], n_perm=n_perm,
             p_upper=float(cs.stats.perm_p(null, obs, "greater", tie=0.0)) if n_perm else np.nan)
    v["H81_1"] = "PASS" if (v["D"] > 0 and v["p_upper"] < 0.05) else "FAIL"
    return e1, v


def a81_sets(g: pd.DataFrame) -> pd.DataFrame:
    """A81 sets from the selection-half table: rhythmic (p_fast < 0.01), eligible nt >= 4 and R >= 0.7, lock tertiles."""
    r = g[(g.p_fast < 0.01) & ~g.gene.isin(G.EXCLUDE)]
    agg = r.groupby("gene").agg(n_rhythmic=("tissue", "size"), amp=("amp_fast", "median"), lock=("lock", "median"),
                                R=("acro_h", G.circ_R),
                                acro_h=("acro_h", lambda x: np.angle(np.exp(1j * W * x).mean()) / W % 24))
    el = agg[(agg.n_rhythmic >= 4) & (agg.R >= 0.7)]
    if not len(el):
        return el.assign(set=pd.Series(dtype=str))
    lo, hi = el.lock.quantile([1 / 3, 2 / 3])
    return pd.concat([el[el.lock <= lo].nlargest(G.N_SET, "amp").assign(set="cue"),
                      el[el.lock >= hi].nlargest(G.N_SET, "amp").assign(set="locked")])


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm = 50 if run["quick"] else 2000
    o = G.obs_table()
    sub = G.pick_donors(o, args.subsample, args.seed)
    H = G.halves(o, sub)
    print({k: {h: len(s) for h, s in v.items()} for k, v in H.items()}, flush=True)

    # ---- A81: selection halves only
    f = out / "a81_selection_genes.csv.gz"
    if f.exists():
        g81 = pd.read_csv(f)
    else:
        g81 = G.selection(G.A81_TISSUES, H["fast"]["sel"], H["vent"]["sel"], H["vent"]["sel"], jobs=args.jobs)
        g81.to_csv(f, index=False)
    sets = a81_sets(g81)
    sets.to_csv(out / "a81_sets.csv")
    print(f"A81: {sets.set.value_counts().to_dict() if len(sets) else 0} genes in the sets", flush=True)
    if len(sets) and sets.set.value_counts().reindex(["cue", "locked"]).fillna(0).min() >= 3:
        wfs = {"core": cs.gtex(), **{k: G.waveform_of(x, k) for k, x in sets.groupby("set")}}
        e1, v = evaluate(eval_tables(G.A81_TISSUES, H, wfs), wfs, n_perm)
        e1.to_csv(out / "a81_e1.csv", index=False)
        pd.Series(v).to_csv(out / "a81_h1_verdict.csv", header=False)
        print(pd.Series(v).to_string(), flush=True)

    # ---- A81b: all fast deaths + ventilated selection half; 26-tissue table also for A89-A91, A100-A102b
    g26 = G.selection_tables(out, args.subsample, args.seed, args.jobs)
    g81b = g26[g26.tissue.isin(G.A81B_TISSUES)]
    agg = G.per_gene(g81b)
    fc = G.frame_check(agg)
    pd.Series(fc).to_csv(out / "a81b_frame_check.csv", header=False)
    print("frame check", fc, flush=True)
    sets_b, med, n_el = G.cue_sets(agg)
    print(f"A81b: eligible {n_el}, median lock {med:.3f}; {sets_b.set.value_counts().to_dict()}", flush=True)
    sets_b.to_csv(out / "a81b_sets.csv")
    ref = sets_b[["set", "nt", "fisher", "Rw", "lock", "amp", "acro_h"]]
    ref.rename_axis("gene").to_csv(out / "a81b_cue_reference.csv")
    (out / "a81b_cue_reference.README.txt").write_text(README.format(
        frame=f"circular correlation {fc['circ_corr']:.2f} over {fc['genes']} genes (guard 0.8: failed, overridden post hoc "
              f"in the legacy run), mean shift {fc['mean_shift_h']:.2f} h"))
    if len(sets_b) and sets_b.set.value_counts().min() >= 15:
        wfs = {"core": cs.gtex(), **{k: G.waveform_of(x, k) for k, x in sets_b.groupby("set")}}
        Hb = {"fast": {"eval": H["fast"]["eval"]}, "vent": {"eval": H["vent"]["eval"]}}
        e1, v = evaluate(eval_tables(G.A81B_TISSUES, Hb, wfs), wfs, 0)
        e1.to_csv(out / "a81b_e1.csv", index=False)
        print(e1.round(2).to_string(index=False), flush=True)
    else:
        print("A81b infeasible: a set has fewer than 15 genes", flush=True)


if __name__ == "__main__":
    main()
