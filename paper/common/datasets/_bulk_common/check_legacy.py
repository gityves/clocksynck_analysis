"""Equivalence checks of the bulk / paired builders against the legacy in-memory constructions (run from anywhere:
`python -I check_legacy.py [name ...]`). Prints one line per cohort / dataset: PASS when rows, clock genes, depth and
(where the legacy has them) expression-matched means agree to numerical precision. Needs the legacy `sccc` install."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # paper/
from common.paths import ARCHIVE, DATASETS  # noqa: E402
LEG = ARCHIVE / "legacy" / "desync" / "validation"
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
sys.path.insert(0, str(LEG / "paired_damping" / "scripts"))
sys.path.insert(0, str(HERE))
import bulk_tables as T  # noqa: E402

GENES = T.GENES
RES = []


def builder(name):
    f = DATASETS / name / "build.py"
    spec = importlib.util.spec_from_file_location(f"b_{name}", f)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def report(name, ok, note=""):
    RES.append((name, "PASS" if ok else "FAIL", note))
    print(f"{RES[-1][1]}  {name}  {note}", flush=True)


def same_pairs(samples, cohort, les, non, counts, tag):
    """samples (mine) vs a legacy (lesional, non-lesional) pair of frames indexed by patient."""
    g = samples[samples.cohort == cohort]
    ml = g[g.condition == "lesional"].set_index("unit")
    mn = g[g.condition == "nonlesional"].set_index("unit")
    pts = sorted(set(ml.index) & set(mn.index))
    ok = list(pts) == list(les.index) == list(non.index)
    if ok:
        for mine, leg in ((ml.loc[pts], les), (mn.loc[pts], non)):
            a, b = mine[GENES].to_numpy(float), leg.rename(columns={"BMAL1": "ARNTL"})[GENES].to_numpy(float)
            ok &= np.allclose(a, b, equal_nan=True)
            if counts:
                ok &= np.allclose(mine.depth.to_numpy(float), leg.depth.to_numpy(float))
        ok &= (mine.value_type.iloc[0] == "counts") == counts
    report(f"{tag} {cohort}", bool(ok), f"{len(pts)} pairs (legacy {len(les)})")


def check_means(mine: dict, legacy_fn, items, tag, key=lambda s: s):
    for item in items:
        gse, les, non, counts = item
        k = key(gse)
        a = mine[k]
        b = legacy_fn(gse, non, counts)
        b = b.rename(index={"BMAL1": "ARNTL"})
        a = a.rename(index={"BMAL1": "ARNTL"})
        ok = a.index.equals(b.index) or set(a.index) == set(b.index)
        ok = ok and np.allclose(a.reindex(b.index).to_numpy(float), b.to_numpy(float), equal_nan=True)
        report(f"{tag} means {k}", bool(ok), f"{len(b)} genes")


def skin(which):
    import a53_skin as A53
    import a54_skin as S
    import a55_ibd as A55
    import a59_other as A59
    if which == "skin_paired_discovery":
        m = builder(which)
        s, _ = m.build()
        for stem, dis, les, non in A53.series():
            same_pairs(s, f"{stem}|{dis}", les, non, les.value_type.iloc[0] == "counts", "a53")
        disc = m.psoriasis_discovery(s)
        items = S.load("discovery")
        for gse, les, non, counts in items:
            same_pairs(disc, f"{gse}|psoriasis", les, non, counts, "a54 discovery")
        report("a54 discovery series", [i[0] for i in items] == sorted(m.DISCOVERY_PSORIASIS) or
               set(i[0] for i in items) == set(m.DISCOVERY_PSORIASIS), f"{len(items)} series")
        gm = m.gene_means()
        check_means(gm, S.match_means, items, "a54 discovery", key=lambda g: f"{g}|psoriasis")
    elif which == "skin_paired_heldout":
        m = builder(which)
        s, _ = m.build()
        items = S.load("heldout")
        for gse, les, non, counts in items:
            same_pairs(s, f"{gse}|psoriasis", les, non, counts, "a54 heldout")
        check_means(m.gene_means(), S.match_means, items, "a54 heldout", key=lambda g: f"{g}|psoriasis")
    elif which == "skin_paired_other":
        m = builder(which)
        s, _ = m.build()
        items = A59.load()
        for name, les, non, counts in items:
            gse, dis = name.split("_", 1)
            same_pairs(s, f"{gse}|{dis}", les, non, counts, "a59")
        lin = A59.load(src=A59.SRC.parent / "processed_other_linear")
        for name, les, non, counts in lin:
            gse, dis = name.split("_", 1)
            same_pairs(s, f"{gse}|{dis}", les, non, counts, "a59 linear")
        check_means(m.gene_means(), S.match_means, items + lin, "a59", key=lambda n: n.replace("_", "|", 1))
    elif which == "gut_bulk_ibd":
        m = builder(which)
        s, _ = m.build()
        items = A55.load()
        for gse, les, non, counts in items:
            same_pairs(s, gse, les, non, counts, "a55")
        for dis in ("CD", "UC"):
            for gse, les, non, counts in A55.load(dis):
                g = s[(s.cohort == gse) & (s.disease == dis)]
                a = set(g[g.condition == "lesional"].unit) & set(g[g.condition == "nonlesional"].unit)
                report(f"a55 {dis} {gse} pairs", len(a) == len(les), f"{len(a)} vs {len(les)}")
        check_means(m.gene_means(), A55.match_means, items, "a55")
    elif which == "skin_paired_treatment":
        import a56_treatment as A56
        m = builder(which)
        s, _ = m.build()
        for f in ("GSE117468.parquet", "GSE201827.parquet", "GSE201827_w52.parquet"):
            t, counts, d = A56.visits(f, None)
            g = s[s.cohort == f"{f[:-8]}|psoriasis"]
            ok = True
            for (visit, status), leg in t.items():
                x = g[(g.visit == visit) & (g.condition == status)].set_index("unit").loc[leg.index]
                ok &= np.allclose(x[GENES].to_numpy(float), leg.rename(columns={"BMAL1": "ARNTL"})[GENES].to_numpy(float))
            report(f"a56 {f}", bool(ok), f"{len(leg)} patients with all four samples")
    elif which == "skin_layers":
        import a62_layers as A62
        m = builder(which)
        s, _ = m.build()
        for layer in ("epidermis", "dermis"):
            x = m.layer_pairs(s, layer)
            for gse in A62.SERIES:
                r = A62.load(gse, layer)
                if r is None:
                    continue
                les, non = r
                g = x[x.cohort == gse]
                ml = g[g.condition == "lesional"].set_index("unit").loc[les.index]
                mn = g[g.condition == "nonlesional"].set_index("unit").loc[non.index]
                ok = (np.allclose(ml[GENES].to_numpy(float), les[GENES].to_numpy(float)) and
                      np.allclose(mn[GENES].to_numpy(float), non[GENES].to_numpy(float)) and
                      np.allclose(ml.depth, les.depth) and np.allclose(mn.depth, non.depth))
                report(f"a62 {layer} {gse}", bool(ok), f"{len(les)} pairs")


def tcga():
    import a52_tcga as C
    m = builder("tcga_paired")
    for c in C.CANCERS:
        ct, cn = C.pairs(c)
        s, _ = m.build([c])
        t = s[s.condition == "tumour"].set_index("unit").loc[ct.index]
        n = s[s.condition == "normal"].set_index("unit").loc[cn.index]
        ok = (np.allclose(t[GENES].to_numpy(float), ct[GENES].to_numpy(float)) and
              np.allclose(n[GENES].to_numpy(float), cn[GENES].to_numpy(float)) and
              np.allclose(t.depth, ct.sum(1)) and np.allclose(n.depth, cn.sum(1)))
        lm = np.log(cn.div(cn.sum(1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
        gm = m.gene_means([c])[c]
        ok &= np.allclose(gm.reindex(lm.index).to_numpy(float), lm.to_numpy(float))
        report(f"a52 {c}", bool(ok), f"{len(ct)} pairs")


def timed_skin():
    sys.path.insert(0, str(LEG / "skin_paired" / "scripts"))
    import skin_paired as SP
    from common.paths import DATA
    m = builder("skin_timed_gse205155")
    s, mm, sym = SP.load(DATA / "skin_paired" / "GSE205155")
    leg = SP.table(s, mm, sym)
    pb, donors = m.build()
    a = pb.set_index(["donor", "compartment"]).sort_index()
    b = leg.set_index(["resection", "compartment"]).sort_index()
    ok = a.index.equals(b.index) and np.allclose(a[GENES].to_numpy(float), b[GENES].to_numpy(float)) and \
        np.allclose(a.depth, b.depth) and (a.stratum.to_numpy() == b.stratum.to_numpy()).all()
    ok &= np.allclose(donors.clock_h.reindex(b.index.get_level_values(0)).to_numpy(float), b.clock_h.to_numpy(float))
    report("a37 GSE205155", bool(ok), f"{len(pb)} rows")


def mouse_ageing():
    names = ("mouse_ageing_wolff", "mouse_ageing_deota")
    built = {n: (builder(n), builder(n).build()[0]) for n in names}  # before the legacy `common` shadows ours
    for k in [k for k in sys.modules if k == "common" or k.startswith("common.")]:
        del sys.modules[k]
    sys.path.insert(0, str(LEG / "mouse_ageing" / "scripts"))
    import a48_wolff as W
    import a49_deota as D
    for name, mod in zip(names, (W, D)):
        m, pb = built[name]
        leg = mod.load()
        a = pb.set_index(["donor", "compartment"]).sort_index()
        b = leg.set_index(["resection", "compartment"]).sort_index()
        ok = a.index.equals(b.index) and np.allclose(a[GENES].to_numpy(float), b[GENES].to_numpy(float)) and \
            np.allclose(a.depth, b.depth)
        report(name, bool(ok), f"{len(pb)} rows, {int(m.tissue_qc(pb).kept.sum())} tissues pass QC")
    for k in [k for k in sys.modules if k == "common" or k.startswith("common.")]:
        del sys.modules[k]
    sys.path.remove(str(LEG / "mouse_ageing" / "scripts"))


def blood():
    for sub in ("blood_clock", "blood_composition", "gonogo", "glucocorticoid"):
        sys.path.insert(0, str(LEG / sub / "scripts"))
    sys.path.insert(0, str(LEG.parents[1] / "scripts"))  # legacy scripts expect <repo>/desync: add their folders
    import blood_clock as BC
    m = builder("blood_timed_gse39445")
    Y, M, ob, wf, kept, nprobe = BC.load()
    s, _ = m.build()
    ok = np.allclose(s[GENES].to_numpy(float) * np.log(2), Y) and np.allclose(s.clock_h.to_numpy(float), ob.clock.to_numpy(float)) \
        and (s.block.to_numpy() == ob.block.to_numpy()).all()
    report("a35 GSE39445", bool(ok), f"{len(s)} samples")


def breast():
    sys.path.insert(0, str(LEG / "scripts"))
    sys.path.insert(0, str(LEG.parent / "common"))
    import v25_5_breast_timed as V
    m = builder("breast_timed_gse233242")
    for typo in (False, True):
        for n in (30e6, 10e6):
            key = V.samples(typo)
            leg = V.pseudobulk(key, n)
            pb, _ = m.build(n_reads=n, join_typo=typo)
            a = pb.set_index(["donor", "compartment"]).sort_index()
            b = leg.set_index(["resection", "compartment"]).sort_index()
            ok = a.index.equals(b.index) and np.allclose(a[GENES].to_numpy(float), b[GENES].to_numpy(float)) and \
                np.allclose(a.depth, b.depth)
            report(f"v25.5 typo={typo} N={n:.0e}", bool(ok), f"{len(pb)} rows")


def baboon_sheep():
    for name, f, key in (("baboon", LEG / "baboon" / "output" / "clock_counts.csv", "resection"),
                         ("sheep", LEG / "farm_gtex" / "output" / "sheep_clock.csv", "resection")):
        m = builder(name)
        pb, _ = m.build()
        cc = pd.read_csv(f)
        cc[key] = cc[key].astype(str)
        if name == "baboon":
            kept = set(pd.read_csv(LEG / "baboon" / "output" / "tissue_qc.csv").query("kept").compartment)
            cc = cc[cc.compartment.isin(kept)]
        a = pb.set_index(["donor", "compartment"]).sort_index()
        b = cc.set_index([key, "compartment"]).sort_index()
        ok = a.index.equals(b.index) and np.allclose(a[GENES].to_numpy(float), b[GENES].to_numpy(float)) and \
            np.allclose(a.depth, b.depth)
        report(f"{name} vs legacy output csv", bool(ok), f"{len(pb)} rows")


def mouse_skin():
    from common.paths import DATA
    m = builder("mouse_skin_timed")
    s, _ = m.build()
    ok = True
    for g in m.SERIES:
        d = pd.read_parquet(DATA / "mouse_skin_timed" / f"timed_{g}.parquet").reset_index(drop=True)
        x = s[s.cohort == g]
        ok &= len(d) == len(x) and np.allclose(d[GENES].to_numpy(float), x[GENES].to_numpy(float), equal_nan=True) and \
            np.allclose(d.zt.astype(float) % 24, x.zt)
    report("a72 mouse skin parquets", bool(ok), f"{len(s)} samples")


STEPS = {"skin": ["skin_paired_discovery", "skin_paired_heldout", "skin_paired_other", "gut_bulk_ibd",
                  "skin_paired_treatment", "skin_layers"]}

if __name__ == "__main__":
    want = sys.argv[1:] or ["skin", "tcga", "timed_skin", "mouse_ageing", "blood", "breast", "baboon_sheep", "mouse_skin"]
    for w in want:
        if w == "skin":
            for n in STEPS["skin"]:
                skin(n)
        elif w in STEPS["skin"]:
            skin(w)
        else:
            globals()[w]()
    print("\n".join(f"{b}\t{a}" for a, b, _ in RES if b == "FAIL") or "all checks PASS")
