"""Cardoso-Moreira 2019: which individual is each matrix column? Prep ONLY, no clock statistic; clock genes are never read.

The expression matrix names its columns `<Organ>.<stage>.<k>` and carries no library name. Supplementary Table 2 of the paper
(`raw/41586_2019_1338_MOESM2_ESM.xlsx`, sheet "SI Table 2") lists every library with the same stage labels, the sex, whether the library was
used (297 of 313, the number of matrix columns) and, for prenatal libraries, the HDBR individual. A column's candidates are therefore the used
libraries of its organ, stage and sex (sex of a column: Y-gene expression, `samples.csv`).

link = exact    one candidate.
link = marker   several candidates, resolved with three deletion-like on/off genes (GSTM1, NPIPB15, COX20P1): a column is given to a library
                when every assignment of its group that contradicts no marker agrees on that library. The marker state of an individual
                comes from its already linked columns in other organs; the rule is iterated to convergence.
link = panel    still several candidates, resolved with a panel of genotype-like genes: autosomal non-clock genes whose expression
                (log2 RPKM, local developmental trend removed within the organ) is concordant between the organs of one individual
                (r >= 0.8 over the same-individual pairs of exact and marker columns; cis-eQTL and copy-number genes such as CCZ1, USP6,
                NPIPA5, GSTT2, PILRB). A column is given to a library when the best assignment of its group beats every assignment that
                places the column elsewhere by >= MARGIN in profile correlation with the individual's linked columns of other organs.
                Columns are only added, never revised. Error rate estimated by leaving each individual out of the panel selection.
link = pair     as panel, for two groups of different organs taken together: the assignment of both groups is scored by the profile
                correlation of the two columns given to each shared individual (plus the evidence from linked columns), same margin.
                This reaches individuals none of whose columns is linked yet.
link = open     not resolved; `library` and `individual` hold the candidates in table order (a guess within the same organ, stage and sex).

The markers are genotypes, not states (on/off calls, or strongly heritable levels): they are validated on the same-individual pairs of
exact columns (tables printed below). An error within a group exchanges two individuals of one organ, stage and sex, which can only weaken
an across-organ statistic whose null permutes individuals within age strata. Clock genes and sex chromosomes are excluded from the panel.

    python -I link_cardoso_individuals.py
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd

DAT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/natexp/cardoso_moreira2019")
OUT = Path(__file__).resolve().parents[1] / "output"
MARKERS = {"ENSG00000134184": "GSTM1", "ENSG00000196436": "NPIPB15", "ENSG00000213025": "COX20P1"}
ORGAN = {"Forebrain/Cerebrum": "Brain", "Hindbrain/Cerebellum": "Cerebellum"}
REF = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/reference/gene_positions_gencode46.tsv")
PANEL_R, MARGIN, WINDOW = 0.8, 0.3, 5   # panel concordance; decision margin; half-width (samples in age order) of the local trend
MIN_GAP = 1.5      # log2 units between the off and on groups of an organ for a marker to be called there
MIN_FRAC = 0.1     # each group holds at least this fraction of the organ's samples
AGE_BINS = [(4, 8), (9, 12), (13, 20)]


def library_table() -> pd.DataFrame:
    t = pd.read_excel(DAT / "raw/41586_2019_1338_MOESM2_ESM.xlsx", sheet_name="SI Table 2", header=2)
    t = t[t["Species (strain)"].astype(str).str.startswith("Human")]
    t = t.rename(columns={"library ID": "library", "Used library?": "used", "Organ": "organ", "Developmental stage": "stage", "Sex": "sex",
                          "Stage (detail)": "stage_detail", "Total number reads": "reads", "Total number aligned reads": "aligned",
                          "RQN": "rqn", "RIN": "rin", "HDBR ID": "individual"})
    t = t[["library", "used", "organ", "stage", "sex", "stage_detail", "reads", "aligned", "rqn", "rin", "individual"]].copy()
    t["library"] = t.library.astype(str)
    t["used"] = t.used.astype(str).str.strip().str.capitalize()
    t["sex"] = t.sex.str.lower()
    t["organ_matrix"] = t.organ.replace(ORGAN)
    t["individual"] = t.individual.map(lambda v: "" if pd.isna(v) else str(int(v)))
    return t.reset_index(drop=True)


def marker_calls(m: pd.DataFrame) -> pd.DataFrame:
    """On/off call per marker and column (NaN where the organ shows no clear gap), thresholds set per organ at the largest gap."""
    x = pd.read_parquet(DAT / "expr_rpkm.parquet").loc[list(MARKERS)]
    lx = np.log2(x + 0.1)
    call = pd.DataFrame(np.nan, index=MARKERS.values(), columns=m.index)
    for o, cols in m.groupby("organ_matrix").groups.items():
        lo = max(2, int(np.ceil(MIN_FRAC * len(cols))))
        for gid, sym in MARKERS.items():
            v = lx.loc[gid, cols]
            s = np.sort(v.to_numpy())
            gaps = np.diff(s)[lo - 1:len(s) - lo]
            j = gaps.argmax() + lo - 1
            if gaps.max() >= MIN_GAP:
                call.loc[sym, cols] = (v > (s[j] + s[j + 1]) / 2).astype(float)
    return call


def candidates(m: pd.DataFrame, t: pd.DataFrame) -> dict[tuple, tuple[list, list]]:
    """Groups of columns and their candidate libraries (same organ, stage and sex; organ and stage only where the sexes disagree)."""
    u = t[t.used == "Yes"]
    groups = {}
    for (o, s), gg in m.groupby(["organ_matrix", "stage_table"]):
        tt = u[(u.organ_matrix == o) & (u.stage == s)]
        assert len(tt) == len(gg), (o, s, len(tt), len(gg))
        if sorted(gg.sex) == sorted(tt.sex):
            for sx in gg.sex.unique():
                groups[o, s, sx] = (gg.index[gg.sex == sx].tolist(), tt.library[tt.sex == sx].tolist())
        else:
            groups[o, s, "any"] = (gg.index.tolist(), tt.library.tolist())
    return groups


def resolve(groups, call: pd.DataFrame, ind: pd.Series) -> pd.DataFrame:
    """Iterative marker resolution; returns one row per column (library, individual, link, n_candidates)."""
    rows = {c: dict(library=libs[0], link="exact", n_candidates=1) for cols, libs in groups.values() if len(cols) == 1 for c in cols}
    while True:
        linked = pd.Series({c: r["library"] for c, r in rows.items()})
        who = ind.reindex(linked.values).to_numpy()
        state = call[linked.index].T.assign(individual=who)
        state = state[state.individual != ""].groupby("individual").agg(lambda v: v.dropna().iloc[0] if v.notna().any() and v.dropna().nunique() == 1 else np.nan)
        new = 0
        for cols, libs in groups.values():
            todo = [c for c in cols if c not in rows]
            if not todo:
                continue
            free = [l for l in libs if l not in {rows[c]["library"] for c in cols if c in rows}]
            ok = []
            for perm in itertools.permutations(free):
                bad = 0
                for c, l in zip(todo, perm):
                    i = ind[l]
                    if i in state.index:
                        a, b = call[c], state.loc[i]
                        bad += int(((a != b) & a.notna() & b.notna()).sum())
                if bad == 0:
                    ok.append(perm)
            for k, c in enumerate(todo):
                opts = {p[k] for p in ok}
                if len(opts) == 1 and len(free) > 1:
                    rows[c] = dict(library=opts.pop(), link="marker", n_candidates=len(libs))
                    new += 1
                elif len(free) == 1:          # the last library of a group whose other columns were resolved
                    rows[c] = dict(library=free[0], link="marker", n_candidates=len(libs))
                    new += 1
        if not new:
            break
    for cols, libs in groups.values():
        todo = [c for c in cols if c not in rows]
        free = [l for l in libs if l not in {rows[c]["library"] for c in cols if c in rows}]
        for c, l in zip(todo, free):
            rows[c] = dict(library=l, link="open", n_candidates=len(libs))
    return pd.DataFrame(rows).T


def profiles(m: pd.DataFrame) -> pd.DataFrame:
    """Gene x column z-scores of log2(RPKM + 0.1) around the local developmental trend of the organ (autosomal non-clock genes)."""
    from clocksync.genes import CLOCK_GENES
    x = pd.read_parquet(DAT / "expr_rpkm.parquet")
    sym = pd.read_csv(DAT / "genes.csv").set_index("gene_id").symbol
    pos = pd.read_csv(REF, sep="\t", header=None, names=["id", "symbol", "chrom", "start", "end", "type"]).set_index("id")
    keep = x.index.intersection(pos.index[~pos.chrom.isin(["chrX", "chrY", "chrM"])])
    keep = keep.difference(sym.index[sym.isin(list(CLOCK_GENES) + ["ARNTL", "BMAL1"])])
    x = x.loc[keep]
    lx = np.log2(x[(x > 1).mean(axis=1) > 0.2] + 0.1)
    order = m["sample"].str.rsplit(".", n=1).str[1].astype(int) if "sample" in m else m.index.str.rsplit(".", n=1).str[1].astype(int)
    k = pd.Series(np.asarray(order), index=m.index)
    z = []
    for _, cols in m.groupby("organ_matrix").groups.items():
        cols = k[cols].sort_values().index
        v, n = lx[cols].to_numpy(), len(cols)
        r = v - np.stack([v[:, max(0, i - WINDOW):min(n, i + WINDOW + 1)].mean(1) for i in range(n)], 1)
        z.append(pd.DataFrame(r / np.maximum(r.std(1, keepdims=True), 1e-6), index=lx.index, columns=cols))
    return pd.concat(z, axis=1)[m.index]


def concordance(Z: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    a, b = Z[:, pairs[:, 0]], Z[:, pairs[:, 1]]
    return (a * b).mean(1) / np.sqrt((a ** 2).mean(1) * (b ** 2).mean(1))


def same_individual_pairs(r: pd.DataFrame, col_pos: dict) -> tuple[np.ndarray, np.ndarray]:
    e = r[r.individual != ""]
    pr = [(col_pos[i], col_pos[j], h) for h, g in e.groupby("individual") for i, j in itertools.combinations(g.index, 2)]
    return np.array([q[:2] for q in pr]), np.array([q[2] for q in pr])


def panel_resolve(groups, r: pd.DataFrame, Zp: pd.DataFrame, ind: pd.Series, m: pd.DataFrame) -> pd.DataFrame:
    """Add columns whose group assignment is decided by the panel with a margin >= MARGIN; iterated, additions only."""
    r = r.copy()
    Zn = (Zp - Zp.mean()) / Zp.std()
    while True:
        linked = r[(r.link != "open") & (r.individual != "")]
        new = 0
        for cols, libs in groups.values():
            todo = [c for c in cols if r.link[c] == "open"]
            if not todo:
                continue
            free = [l for l in libs if l not in set(r.library[[c for c in cols if c not in todo]])]

            def score(c, l):
                o = linked.index[(linked.individual == ind[l]) & (linked.organ_matrix != m.organ_matrix[c])]
                return np.nan if ind[l] == "" or not len(o) else float((Zn[c] * Zn[o].mean(axis=1)).mean() / Zn[o].mean(axis=1).std())

            S = np.array([[score(c, l) for l in free] for c in todo])
            S0 = np.nan_to_num(S)
            tot = sorted(((sum(S0[i, q[i]] for i in range(len(todo))), q) for q in itertools.permutations(range(len(free)))), reverse=True)
            best = tot[0][1]
            for i, c in enumerate(todo):
                alt = [t for t, q in tot if q[i] != best[i]]
                if alt and not np.isnan(S[i, best[i]]) and tot[0][0] - alt[0] >= MARGIN:
                    r.loc[c, ["library", "individual", "link", "margin"]] = [free[best[i]], ind[free[best[i]]], "panel", tot[0][0] - alt[0]]
                    new += 1
            todo = [c for c in cols if r.link[c] == "open"]          # keep the remaining guesses a bijection onto the remaining libraries
            free = [l for l in libs if l not in set(r.library[[c for c in cols if c not in todo]])]
            for c, l in zip(todo, free):
                r.loc[c, ["library", "individual"]] = [l, ind[l]]
        if not new:
            return r


def pair_resolve(groups, r: pd.DataFrame, Zp: pd.DataFrame, ind: pd.Series, m: pd.DataFrame) -> pd.DataFrame:
    """Joint assignment of two open groups of different organs; additions only, iterated together with `panel_resolve`."""
    Zn = (Zp - Zp.mean()) / Zp.std()
    C = np.corrcoef(Zn.to_numpy().T)
    pos = {c: i for i, c in enumerate(Zn.columns)}
    while True:
        r = panel_resolve(groups, r, Zp, ind, m)
        linked = r[(r.link != "open") & (r.individual != "")]

        def state(cols, libs):
            todo = [c for c in cols if r.link[c] == "open"]
            free = [l for l in libs if l not in set(r.library[[c for c in cols if c not in todo]])]
            lin = np.zeros((len(todo), len(free)))
            for i, c in enumerate(todo):
                for j, l in enumerate(free):
                    o = linked.index[(linked.individual == ind[l]) & (linked.organ_matrix != m.organ_matrix[c])]
                    if ind[l] != "" and len(o):
                        lin[i, j] = np.corrcoef(Zn[c], Zn[o].mean(axis=1))[0, 1]
            return todo, free, lin

        new, keys = 0, [k for k, (cols, _) in groups.items() if any(r.link[c] == "open" for c in cols)]
        for ka, kb in itertools.combinations(keys, 2):
            if ka[0] == kb[0]:
                continue
            ta, fa, la = state(*groups[ka])
            tb, fb, lb = state(*groups[kb])
            shared = {ind[l] for l in fa} & {ind[l] for l in fb} - {""}
            if not ta or not tb or not shared:
                continue
            tot = []
            for qa in itertools.permutations(range(len(fa))):
                wa = {ind[fa[j]]: ta[i] for i, j in enumerate(qa)}
                sa = sum(la[i, j] for i, j in enumerate(qa))
                for qb in itertools.permutations(range(len(fb))):
                    wb = {ind[fb[j]]: tb[i] for i, j in enumerate(qb)}
                    tot.append((sa + sum(lb[i, j] for i, j in enumerate(qb)) + sum(C[pos[wa[h]], pos[wb[h]]] for h in shared), qa, qb))
            tot.sort(key=lambda t: -t[0])
            best = tot[0]
            for side, (todo, free) in enumerate(((ta, fa), (tb, fb))):
                for i, c in enumerate(todo):
                    lib = free[best[1 + side][i]]
                    alt = [t[0] for t in tot if t[1 + side][i] != best[1 + side][i]]
                    if alt and ind[lib] in shared and best[0] - alt[0] >= MARGIN and r.link[c] == "open":
                        r.loc[c, ["library", "individual", "link", "margin"]] = [lib, ind[lib], "pair", best[0] - alt[0]]
                        new += 1
            for cols, libs in (groups[ka], groups[kb]):          # remaining guesses stay a bijection
                todo = [c for c in cols if r.link[c] == "open"]
                free = [l for l in libs if l not in set(r.library[[c for c in cols if c not in todo]])]
                for c, l in zip(todo, free):
                    r.loc[c, ["library", "individual"]] = [l, ind[l]]
        if not new:
            return r


def pair_error(r: pd.DataFrame, Z: pd.DataFrame, m: pd.DataFrame, rng) -> pd.DataFrame:
    """Leave-two-individuals-out for the pair rule: individuals A and B of a close age with columns in the same two organs;
    d = (corr of A's columns + corr of B's columns) - (corr of the two crossed pairs), panel chosen without A and B."""
    col_pos = {c: i for i, c in enumerate(Z.columns)}
    e = r[(r.link.isin(["exact", "marker"])) & (r.individual != "") & r.prenatal]
    pairs, who = same_individual_pairs(e, col_pos)
    Zv, out = Z.to_numpy(), []
    tab = e.reset_index().pivot_table(index="individual", columns="organ_matrix", values="sample", aggfunc="first")
    age = e.groupby("individual").age_wpc.first()
    quads = [(a, b, o1, o2) for a, b in itertools.combinations(tab.index, 2) if abs(age[a] - age[b]) <= 2
             for o1, o2 in itertools.combinations(tab.columns, 2) if tab.loc[[a, b], [o1, o2]].notna().all().all()]
    for k in rng.permutation(len(quads))[:400]:
        a, b, o1, o2 = quads[k]
        top = concordance(Zv, pairs[(who != a) & (who != b)]) >= PANEL_R
        cols = [tab.loc[a, o1], tab.loc[a, o2], tab.loc[b, o1], tab.loc[b, o2]]
        c = np.corrcoef(Z.loc[top, cols].to_numpy().T)
        out.append(dict(a=a, b=b, organs=f"{o1}+{o2}", n_panel=int(top.sum()), d=c[0, 1] + c[2, 3] - c[0, 3] - c[2, 1]))
    return pd.DataFrame(out)


def panel_error(r: pd.DataFrame, Z: pd.DataFrame, m: pd.DataFrame) -> pd.DataFrame:
    """Leave-one-individual-out: a linked column against its own individual and against each other individual of a close age (within 2 wpc,
    other organs), with the panel chosen without the left-out individual. d = correlation with own minus with the other."""
    col_pos = {c: i for i, c in enumerate(Z.columns)}
    e = r[(r.link != "open") & (r.individual != "") & r.prenatal]
    pairs, who = same_individual_pairs(e, col_pos)
    Zv, out = Z.to_numpy(), []
    for h in np.unique(who):
        top = concordance(Zv, pairs[who != h]) >= PANEL_R
        Zn = Z[top]
        Zn = (Zn - Zn.mean()) / Zn.std()
        mine = e.index[e.individual == h]
        for c in mine:
            own = [o for o in mine if m.organ_matrix[o] != m.organ_matrix[c]]
            if not own:
                continue
            cor = lambda o: float(np.corrcoef(Zn[c], Zn[o].mean(axis=1))[0, 1])
            a = cor(own)
            for h2, g in e[(e.individual != h) & ((e.age_wpc - e.age_wpc[c]).abs() <= 2) & (e.organ_matrix != m.organ_matrix[c])].groupby("individual"):
                out.append(dict(column=c, individual=h, other=h2, n_panel=int(top.sum()), d=a - cor(list(g.index))))
    return pd.DataFrame(out)


def validate(r: pd.DataFrame, call: pd.DataFrame, m: pd.DataFrame, link: list[str]) -> pd.DataFrame:
    """Marker agreement over same-individual column pairs (other organs) against the agreement expected by chance."""
    e = r[r.link.isin(link) & (r.individual != "")]
    pairs = [(i, j) for _, g in e.groupby("individual") for i, j in itertools.combinations(g.index, 2)]
    out = []
    for sym in call.index:
        a = call.loc[sym, [i for i, _ in pairs]].to_numpy()
        b = call.loc[sym, [j for _, j in pairs]].to_numpy()
        k = ~np.isnan(a) & ~np.isnan(b)
        f = call.loc[sym, m.index[m.prenatal]].mean()
        out.append(dict(marker=sym, link="+".join(link), pairs=int(k.sum()), agree=int((a[k] == b[k]).sum()), frac_on=round(f, 3),
                        chance=round(f ** 2 + (1 - f) ** 2, 3)))
    return pd.DataFrame(out)


def main():
    t = library_table()
    t.to_csv(DAT / "library_table.csv", index=False)
    m = pd.read_csv(DAT / "samples.csv").set_index("sample")
    m["stage_table"] = m.stage.replace({"Senior": "senior"})
    print("libraries", len(t), "used", (t.used == "Yes").sum(), "| matrix columns", len(m))
    call = marker_calls(m)
    print("columns with a call per marker:", call.notna().sum(axis=1).to_dict(), "| on fraction", call.mean(axis=1).round(2).to_dict())
    groups = candidates(m, t)
    ind = t.set_index("library").individual
    r = resolve(groups, call, ind)
    r["individual"] = ind.reindex(r.library).to_numpy()
    r = r.join(t.set_index("library")[["reads", "aligned", "rqn", "rin", "stage_detail"]], on="library").reindex(m.index)
    r = r.join(m[["organ", "organ_matrix", "stage", "age_bin", "prenatal", "age_wpc", "sex"]])
    r["margin"] = np.nan

    Z = profiles(m)
    col_pos = {c: i for i, c in enumerate(Z.columns)}
    pairs, _ = same_individual_pairs(r[(r.link != "open") & r.prenatal], col_pos)
    rg = pd.Series(concordance(Z.to_numpy(), pairs), index=Z.index)
    sym = pd.read_csv(DAT / "genes.csv").set_index("gene_id").symbol
    panel = rg[rg >= PANEL_R].sort_values(ascending=False)
    pd.DataFrame({"symbol": sym.reindex(panel.index), "concordance": panel.round(3)}).to_csv(OUT / "link_cardoso_panel.csv")
    print(f"panel: {len(panel)} genes with concordance >= {PANEL_R} over {len(pairs)} same-individual pairs:", ", ".join(sym.reindex(panel.index).fillna("?")))
    err = panel_error(r, Z, m)
    err.to_csv(OUT / "link_cardoso_panel_loo.csv", index=False)
    print(f"leave-one-individual-out, {err.column.nunique()} columns x other individuals of a close age = {len(err)} decisions: "
          f"right with margin {np.mean(err.d >= MARGIN):.3f}, undecided {np.mean(err.d.abs() < MARGIN):.3f}, wrong with margin {np.mean(err.d <= -MARGIN):.4f}; "
          f"panel size {err.n_panel.min()}-{err.n_panel.max()}")
    err2 = pair_error(r, Z, m, np.random.default_rng(0))
    err2.to_csv(OUT / "link_cardoso_pair_loo.csv", index=False)
    print(f"leave-two-individuals-out, pair rule, {len(err2)} decisions: right with margin {np.mean(err2.d >= MARGIN):.3f}, "
          f"undecided {np.mean(err2.d.abs() < MARGIN):.3f}, wrong with margin {np.mean(err2.d <= -MARGIN):.4f}")
    r = pair_resolve(groups, r, Z.loc[panel.index], ind, m)
    r.index.name = "sample"
    r.to_csv(DAT / "individuals.csv")

    val = pd.concat([validate(r, call, m, ["exact"]), validate(r, call, m, ["exact", "marker"]), validate(r, call, m, ["exact", "marker", "panel"]),
                     validate(r, call, m, ["exact", "marker", "panel", "pair"])])
    val.to_csv(OUT / "link_cardoso_marker_validation.csv", index=False)
    print(val.to_string(index=False))
    summ = r.groupby(["prenatal", "link"]).size().unstack(fill_value=0)
    summ.to_csv(OUT / "link_cardoso_summary.csv")
    print(summ.to_string())

    p = r[r.prenatal & (r.link != "open")]
    n_org = p.groupby("individual").organ.nunique()
    rows = []
    for lo, hi in AGE_BINS + [(4, 20)]:
        q = p[p.age_wpc.between(lo, hi)]
        k = q.groupby("individual").organ.nunique()
        row = dict(age_bin=f"{lo}-{hi}", columns=len(q), individuals=k.size, ge2=int((k >= 2).sum()), ge3=int((k >= 3).sum()))
        tab = q.assign(v=1).pivot_table(index="individual", columns="organ", values="v", aggfunc="max", fill_value=0)
        for a, b in itertools.combinations(sorted(tab.columns), 2):
            row[f"{a}+{b}"] = int((tab[a] & tab[b]).sum())
        rows.append(row)
    d = pd.DataFrame(rows).fillna(0)
    d.to_csv(OUT / "design_cardoso_linked.csv", index=False)
    print(d.T.to_string())
    print("linked prenatal individuals by number of organs:", n_org.value_counts().sort_index().to_dict())
    full = r[r.prenatal].groupby("individual").organ.nunique()
    print("all prenatal columns (open ones at their guess): individuals", full.size, "with >= 2 organs", int((full >= 2).sum()))


if __name__ == "__main__":
    main()
