"""A256: does the within-person rotation of the clock state replicate in fresh studies, also when nothing is done?

Registered in paper/docs/registered_experiments.md (A256 and its design addendum). Statistic of A252 (imported):
rho = sum Re(delta conj(c)) / sum |c|^2 with delta = z_k - z_0 per unit and c the predicted chord; null: time labels
permuted within unit; companions: matched and responsive pseudo-clock sets.

Stages
  power  units, null SD and MDE of rho for H256.1 and H256.2 (permuted labels only) -> output/{power,design}.csv
  score  (after the power note is dated) everything registered -> output/score_{tests,cells,descriptive}.csv
  `--placebo`: time labels shuffled once within unit before anything else (code test; output/placebo/).

    python compute_a256.py --stage power|score [--smoke] [--placebo]

Choices not fixed by the registration
- Input: `output/prep/<GSE>_{expr.parquet,sheet.csv}` (`prep_a256.py`); y = ln(count + 1) for counts, deposited log values x ln 2
  for arrays. Classes as in the design addendum (`CLASS`).
- A cell = dataset x arm x later time point. H256.3 counts cells whose implied midpoint (direction of the mean delta) lies
  within 2 h of the protocol's; null: the same count under the label permutations.
- H256.4: D = rho(intervention arms) - rho(nothing done), pooled sums over GSE40589 and GSE59088; null: arm labels permuted
  within person (GSE40589: a person's units exchange arms) or over persons (GSE59088); two-sided on |D - null mean|.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A252_elapsed_time/scripts"))
from clocksync import paired as P  # noqa: E402
from common import cli  # noqa: E402
import compute_a252 as E  # noqa: E402

OUT = HERE.parent / "output"
PREP = HERE.parent / "output" / "prep"
SEED, N_PERM, N_SETS, MDE_K = 256, 10_000, 200, 2.49
WF, W = E.WF, E.W
DSETS = ["GSE40589", "GSE9105", "GSE59088", "GSE220899", "GSE202295", "GSE163356", "GSE126296"]
RECOVERY = {("GSE59088", "endurance"), ("GSE59088", "resistance"), ("GSE202295", "NGT"), ("GSE202295", "T2D")}


class Data:
    def __init__(self, acc, placebo=False):
        self.acc = acc
        x = pd.read_parquet(PREP / f"{acc}_expr.parquet").astype(float)
        s = pd.read_csv(PREP / f"{acc}_sheet.csv", dtype={"person": str, "unit": str}).set_index("sample")
        s = s[s.index.isin(x.index) & (s.elapsed_h <= 12)]
        ok = s.groupby("unit").elapsed_h.transform(lambda v: (v == 0).sum() == 1 and (v > 0).sum() >= 1 and v.nunique() == len(v))
        s = s[ok]
        x = x.loc[s.index].dropna(axis=1)
        self.y = np.log(x + 1.0) if s.kind.astype(str).str.startswith("counts").all() else x * np.log(2.0)
        self.sh = s
        self.cls = np.array(["recovery" if (acc, a) in RECOVERY else c for a, c in zip(s.arm, s.arm_class)])
        self.unit = pd.factorize(s.unit)[0]
        self.el = s.elapsed_h.to_numpy(float).copy()
        if placebo:
            rng = np.random.default_rng([SEED, 9, sum(map(ord, acc))])
            for u in np.unique(self.unit):
                i = np.flatnonzero(self.unit == u)
                self.el[i] = self.el[rng.permutation(i)]
        self.t0 = s.t0_clock_h.to_numpy(float)
        self.wf = WF.subset([g for g in WF.genes if g in self.y.columns])
        self.idx = [np.flatnonzero(self.unit == u) for u in np.unique(self.unit)]

    def parts(self, z, el=None, mask=None):
        """(sum Re(delta conj c), sum |c|^2) over later samples (of the units in `mask`)."""
        el = self.el if el is None else el
        zb = np.empty(len(z), complex)
        for i in self.idx:
            zb[i] = z[i][el[i] == 0][0]
        c = np.exp(1j * W * (self.t0 + el)) - np.exp(1j * W * self.t0)
        m = el > 0 if mask is None else (el > 0) & mask
        return float(np.sum(((z - zb)[m] * np.conj(c[m])).real)), float(np.sum(abs(c[m]) ** 2))

    def perm_el(self, rng):
        el = self.el.copy()
        for i in self.idx:
            el[i] = self.el[rng.permutation(i)]
        return el

    def z(self, wf=None):
        return E.states(self.y, wf or self.wf)

    def pseudo(self, kind, n, tag):
        lm = self.y.mean(0).astype(float)
        return P.pseudo_waveforms(lm[~lm.index.duplicated()], n, np.random.default_rng([SEED, tag, sum(map(ord, self.acc))]), self.wf, kind)

    def cells(self, z, el=None):
        """(arm, elapsed) -> (implied midpoint h, protocol midpoint h, rho, units)."""
        el = self.el if el is None else el
        zb = np.empty(len(z), complex)
        for i in self.idx:
            zb[i] = z[i][el[i] == 0][0]
        d, out = z - zb, {}
        for a in pd.unique(self.sh.arm):
            for e in sorted(set(self.el[self.el > 0])):
                m = (self.sh.arm.to_numpy() == a) & (el == e)
                if m.sum() < 3:
                    continue
                c = np.exp(1j * W * (self.t0[m][0] + e)) - np.exp(1j * W * self.t0[m][0])
                mean = d[m].mean()
                out[(a, e)] = (((np.angle(mean) - np.pi / 2) / W) % 24, (self.t0[m][0] + e / 2) % 24, float(np.sum((d[m] * np.conj(c)).real) / (m.sum() * abs(c) ** 2)), int(m.sum()))
        return out


def rho(parts):
    den = sum(p[1] for p in parts)
    return sum(p[0] for p in parts) / den if den > 0 else np.nan


def within2(cells):
    return int(sum(abs((v[0] - v[1] + 12) % 24 - 12) <= 2 for v in cells.values()))


def p_up(nul, obs):
    return float((1 + np.sum(np.asarray(nul) >= obs)) / (1 + len(nul)))


def run(args):
    n_perm, n_sets = (300, 5) if args.smoke else (N_PERM, N_SETS)
    D = {a: Data(a, args.placebo) for a in DSETS}
    Z = {a: D[a].z() for a in DSETS}
    nd = {a: D[a].cls == "nothing_done" for a in DSETS}
    fam = {"H256.1 nothing done": lambda a, z, el=None: D[a].parts(z, el, nd[a]) if nd[a].any() else (0.0, 0.0),
           "H256.2 all fresh units": lambda a, z, el=None: D[a].parts(z, el)}
    rng = np.random.default_rng([SEED, 1])
    nul = {k: [] for k in fam}
    ncount = []
    for _ in range(n_perm):
        els = {a: D[a].perm_el(rng) for a in DSETS}
        for k, f in fam.items():
            nul[k].append(rho([f(a, Z[a], els[a]) for a in DSETS]))
        ncount.append(sum(within2(D[a].cells(Z[a], els[a])) for a in DSETS))
    design = pd.DataFrame([dict(dataset=a, arm=arm, cls=("recovery" if (a, arm) in RECOVERY else g.arm_class.iloc[0]), units=g.unit.nunique(),
                                later_points=";".join(f"{x:.2f}" for x in sorted(set(g.elapsed_h[g.elapsed_h > 0]))), t0=g.t0_clock_h.iloc[0],
                                clock_genes=len(D[a].wf.genes)) for a in DSETS for arm, g in D[a].sh.groupby("arm")])
    design.to_csv(OUT / "design.csv", index=False)
    power = pd.DataFrame([dict(family=k, units=int(design.units[design.cls == "nothing_done"].sum()) if "nothing" in k else int(design.units.sum()),
                               null_mean=float(np.mean(v)), null_sd=float(np.std(v)), mde=MDE_K * float(np.std(v))) for k, v in nul.items()])
    power["scored"] = [bool(power.units[0] >= 10 and power.mde[0] <= 1.0), bool(power.mde[1] <= 0.7)]
    power.to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 250)
    print(design.to_string(), power.round(4).to_string(), sep="\n", flush=True)
    if args.stage == "power":
        return

    comp = {kind: {a: [D[a].z(w) for w in D[a].pseudo(kind, n_sets, tag)] for a in DSETS} for kind, tag in (("matched", 3), ("responsive", 4))}
    tests = []
    for k, f in fam.items():
        obs = rho([f(a, Z[a]) for a in DSETS])
        row = dict(family=k, rho=obs, null_mean=float(np.mean(nul[k])), null_sd=float(np.std(nul[k])), p=p_up(nul[k], obs), mde=MDE_K * float(np.std(nul[k])))
        for kind in comp:
            v = np.array([rho([f(a, comp[kind][a][j]) for a in DSETS]) for j in range(n_sets)])
            row.update({f"{kind}_median": float(np.median(v)), f"{kind}_q95": float(np.quantile(v, .95)), f"p_{kind}": p_up(v, obs)})
        row["clock_specific"] = bool(obs > row["matched_q95"] and obs > row["responsive_q95"])
        tests.append(row)
    T = pd.DataFrame(tests)
    T["p_holm"] = np.minimum(1, [max(2 * T.p.min(), T.p.max()) if p == T.p.max() and len(set(T.p)) > 1 else 2 * T.p.min() for p in T.p])
    cells = {a: D[a].cells(Z[a]) for a in DSETS}
    n_in, n_cells = sum(within2(c) for c in cells.values()), sum(len(c) for c in cells.values())
    desc = [dict(item="H256.3 cells with implied midpoint within 2 h of the protocol's", value=n_in, cells=n_cells, null_mean=float(np.mean(ncount)), p=p_up(ncount, n_in))]
    C = pd.DataFrame([dict(dataset=a, arm=arm, cls=("recovery" if (a, arm) in RECOVERY else D[a].sh.arm_class[D[a].sh.arm == arm].iloc[0]), elapsed_h=e, units=v[3],
                           implied_midpoint_h=v[0], protocol_midpoint_h=v[1], rho=v[2]) for a in DSETS for (arm, e), v in cells[a].items()])
    # H256.4
    rng = np.random.default_rng([SEED, 6])
    def Dstat(lab):
        pi, pn = [], []
        for a in ("GSE40589", "GSE59088"):
            pi.append(D[a].parts(Z[a], None, lab[a] == "x"))
            pn.append(D[a].parts(Z[a], None, lab[a] == "n"))
        return rho(pi) - rho(pn)
    lab = {a: np.where(D[a].cls == "nothing_done", "n", "x") for a in ("GSE40589", "GSE59088")}
    obs4, n4 = Dstat(lab), []
    for _ in range(2000 if not args.smoke else 100):
        pl = {}
        for a in lab:
            d = D[a]
            ul = np.array([lab[a][i][0] for i in d.idx])
            if a == "GSE40589":
                per = d.sh.person.to_numpy()
                up = np.array([per[i][0] for i in d.idx])
                for q in np.unique(up):
                    j = np.flatnonzero(up == q)
                    ul[j] = ul[rng.permutation(j)]
            else:
                ul = rng.permutation(ul)
            pl[a] = ul[d.unit]
        n4.append(Dstat(pl))
    n4 = np.array(n4)
    desc.append(dict(item="H256.4 rho(intervention or recovery arms) - rho(nothing done), GSE40589 + GSE59088", value=obs4, null_mean=float(n4.mean()), null_sd=float(n4.std()),
                     p=float((1 + np.sum(abs(n4 - n4.mean()) >= abs(obs4 - n4.mean()))) / (1 + len(n4)))))
    for name, mask in (("nothing done + recovery arms (letter of the registration)", {a: np.isin(D[a].cls, ["nothing_done", "recovery"]) for a in DSETS}),
                       ("recovery arms", {a: D[a].cls == "recovery" for a in DSETS}), ("intervention arms", {a: D[a].cls == "intervention" for a in DSETS})):
        desc.append(dict(item=f"rho, {name}", value=rho([D[a].parts(Z[a], None, mask[a]) for a in DSETS if mask[a].any()])))
    wf2 = {a: D[a].wf.drop(["PER1", "NR1D1"]) for a in DSETS}
    Z2 = {a: D[a].z(wf2[a]) for a in DSETS}
    for k, f in fam.items():
        desc.append(dict(item=f"rho without PER1 and NR1D1, {k}", value=rho([f(a, Z2[a]) for a in DSETS])))
    for a in DSETS:
        desc.append(dict(item=f"rho, {a}", value=rho([D[a].parts(Z[a])]), units=len(D[a].idx)))
    T.to_csv(OUT / "score_tests.csv", index=False)
    C.to_csv(OUT / "score_cells.csv", index=False)
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    print(T.round(4).T.to_string(), C.round(3).to_string(), pd.DataFrame(desc).round(4).to_string(), sep="\n", flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--placebo", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args(argv)
    cli.setup(args)
    global OUT
    OUT = OUT / ("placebo" if args.placebo else "smoke" if args.smoke else "")
    OUT.mkdir(parents=True, exist_ok=True)
    run(args)


if __name__ == "__main__":
    main()
