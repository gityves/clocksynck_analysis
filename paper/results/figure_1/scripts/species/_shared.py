"""Helpers shared by the scripts of species. Clock quantities come from `clocksync`; what is local is marked
"local: not in clocksync" and listed in the README."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from common import data  # noqa: E402

W = 2 * np.pi / 24
HOURS = 24 / (2 * np.pi)


def _psi_step(e2, v, psi, floor=1e-8, n_iter=100, tol=1e-10) -> float:
    for _ in range(n_iter):
        w2 = 1.0 / (psi + v) ** 2
        new = max(float(np.sum(w2 * (e2 - v)) / np.sum(w2)), floor)
        if abs(new - psi) <= tol * max(psi, floor):
            return new
        psi = new
    return psi


def timed_fit(beta, V, t_h, tissue, fit_rows=None, n_iter: int = 200, tol: float = 1e-10):
    """local: not in clocksync. Per tissue beta = mu_k + a_k e^{i w t} + e by Poisson-weighted least squares with the
    own-noise variance psi_k by maximum likelihood (the A36 / 01_gtex_organs timed fit). Rows `fit_rows` (default
    all) enter the fit; returns (labels, mu (K x 2), a (complex K), psi, complex residual of every row)."""
    beta, V, t_h = np.asarray(beta, float), np.asarray(V, float), np.asarray(t_h, float)
    fit_rows = np.ones(len(beta), bool) if fit_rows is None else np.asarray(fit_rows)
    lab = np.unique(tissue[fit_rows])
    f = np.column_stack([np.cos(t_h * W), np.sin(t_h * W)])
    mu, a, psi = np.zeros((len(lab), 2)), np.zeros(len(lab), complex), np.zeros(len(lab))
    I2 = np.eye(2)
    J = np.array([[0.0, -1.0], [1.0, 0.0]])
    for kk, k in enumerate(lab):
        i = (tissue == k) & fit_rows
        b, Vi, fi = beta[i], V[i], f[i]
        G = np.zeros((i.sum(), 2, 4))
        G[:, :, :2] = I2
        G[:, :, 2], G[:, :, 3] = fi, fi @ J.T
        v_eig, U = np.linalg.eigh(Vi)
        p = max(float(np.var(b, 0).mean()), 1e-8)
        for _ in range(n_iter):
            Wt = np.linalg.inv(p * I2 + Vi)
            GtW = np.swapaxes(G, 1, 2) @ Wt
            th = np.linalg.solve((GtW @ G).sum(0), np.einsum("nab,nb->a", GtW, b))
            e = b - np.einsum("nab,b->na", G, th)
            new = _psi_step(np.einsum("nai,na,nb,nbi->ni", U, e, e, U).ravel(), v_eig.ravel(), p)
            done = abs(new - p) <= tol * max(p, 1e-8)
            p = new
            if done:
                break
        mu[kk], psi[kk], a[kk] = th[:2], p, th[2] + 1j * th[3]
    k = np.searchsorted(lab, tissue)
    e = (beta[:, 0] + 1j * beta[:, 1]) - (mu[k, 0] + 1j * mu[k, 1]) - a[k] * np.exp(1j * t_h * W)
    return lab, mu, a, psi, e


def known_host(beta, V, f, grp, n_iter: int = 200, tol: float = 1e-10) -> dict:
    """local: not in clocksync. Coupling of one compartment to a host state that is known exactly (here cos / sin of the sampling
    time): the regression of `anchored.on_host` with the A17 small-sample corrections (|z|^2 minus the sampling variance of z,
    psi raised by the degrees of freedom). beta (n, 2), V (n, 2, 2), f (n, 2), grp: intercept group of each row.
    Returns z, z2, psi, V (host variance), coupling, boundary."""
    beta, V, f = np.asarray(beta, float), np.asarray(V, float), np.asarray(f, float)
    gl, gi = np.unique(np.asarray(grp).astype(str), return_inverse=True)
    n = len(beta)
    I2 = np.eye(2)
    J = np.array([[0.0, -1.0], [1.0, 0.0]])
    G = np.zeros((n, 2, 2 * len(gl)))
    for j in range(len(gl)):
        G[gi == j, :, 2 * j:2 * j + 2] = I2
    G = np.concatenate([G, np.stack([f, f @ J.T], -1)], -1)
    v_eig, U = np.linalg.eigh(V)
    p = max(float(np.var(beta, 0).mean()), 1e-8)
    for _ in range(n_iter):
        Wt = np.linalg.inv(p * I2 + V)
        GtW = np.swapaxes(G, 1, 2) @ Wt
        A_ = GtW @ G
        th = np.linalg.solve(A_.sum(0), np.einsum("nab,nb->a", GtW, beta))
        e = beta - np.einsum("nab,b->na", G, th)
        new = _psi_step(np.einsum("nai,na,nb,nbi->ni", U, e, e, U).ravel(), v_eig.ravel(), p)
        done = abs(new - p) <= tol * max(p, 1e-8)
        p = new
        if done:
            break
    Wt = np.linalg.inv(p * I2 + V)
    cov = np.linalg.inv((np.swapaxes(G, 1, 2) @ Wt @ G).sum(0))
    z = th[-2] + 1j * th[-1]
    npar, m = G.shape[2], 2 * n
    vbar = float(np.trace(V, axis1=1, axis2=2).mean() / 2)
    psi = max(p + npar / (m - npar) * (p + vbar), 1e-8) if m > npar else np.nan
    z2 = max(abs(z) ** 2 - float(np.trace(cov[-2:, -2:])), 0.0)
    fc = f.copy()
    for j in range(len(gl)):
        fc[gi == j] -= fc[gi == j].mean(0)
    Vh = float(np.trace(np.cov(fc.T)))
    s = z2 * Vh
    return dict(z=z, z2=z2, psi=psi, V=Vh, coupling=s / (s + 2 * psi), boundary=p <= 1e-7)


def gtex_gains(tissues, wf=None, n_boot: int = 200, subsample=None, seed: int = 0) -> pd.DataFrame:
    """local: not in clocksync. Timed GTEx gain g_t of each tissue (fast deaths with a time of death): cosinor of log(CPM + 0.5) of
    the waveform genes (median CPM >= 5) with the D2 covariates (sex, age bracket, Hardy, ischaemic time and square), regressed on
    the waveform: g = <xc, c> / <xc, xc>; donor bootstrap se. Aggregates only (GTEx controlled access)."""
    from common import gtex

    wf = cs.gtex() if wf is None else wf
    xc = wf.amp * np.exp(1j * wf.psi)
    pb, donors = gtex.load(list(tissues), genes=list(wf.genes), depth="library_size")

    def gain(cpm, t, Z):
        ex = np.median(cpm, 0) >= 5
        c = cs.timed.cosinor(np.log(cpm[:, ex] + 0.5), t, covariates=Z)
        c = c["amplitude"] * np.exp(1j * W * c["acrophase_h"])
        x = xc[ex]
        return complex(np.conj(x) @ c / (np.conj(x) @ x)), int(ex.sum())

    rows = []
    for tis in tissues:
        q = pb[pb.compartment == tis].set_index("donor")
        d = donors.loc[q.index]
        if subsample and subsample < len(q):
            keep = np.random.default_rng(seed).permutation(len(q))[:subsample]
            q, d = q.iloc[keep], d.iloc[keep]
        cpm = q[list(wf.genes)].to_numpy(float) / q.depth.to_numpy(float)[:, None] * 1e6
        t = d.tod.to_numpy(float) % 24
        Z = pd.get_dummies(pd.DataFrame(dict(sex=d.male.astype(str), age=d.age.astype(str), hardy=d.hardy.astype(str))),
                           drop_first=True).to_numpy(float)
        isch = d.ischaemic_h.astype(float).fillna(d.ischaemic_h.median()).to_numpy()
        isch = (isch - isch.mean()) / (isch.std() + 1e-9)
        Z = np.c_[Z, isch, isch ** 2]
        g, ng = gain(cpm, t, Z)
        rng = np.random.default_rng(0)
        bs = []
        for _ in range(n_boot):
            i = rng.choice(len(t), len(t), replace=True)
            try:
                bs.append(gain(cpm[i], t[i], Z[i])[0])
            except np.linalg.LinAlgError:
                continue
        bs = np.array(bs)
        dphi = np.angle(bs * np.conj(g)) / W if len(bs) > 1 else np.array([np.nan, np.nan])
        rows.append(dict(tissue=tis, n=len(t), clock_genes=ng, g_re=g.real, g_im=g.imag, abs_g=abs(g),
                         abs_g_se=np.std(np.abs(bs), ddof=1) if len(bs) > 1 else np.nan, phi_g_h=(np.angle(g) / W) % 24,
                         phi_g_se_h=np.std(dphi, ddof=1) if len(bs) > 1 else np.nan))
    return pd.DataFrame(rows)


def timed_gain(sub: pd.DataFrame, zt, wf) -> tuple:
    """local: not in clocksync. Per-tissue cosinor on ZT of log(CPM + 0.5) of the waveform genes with median CPM >= 5, regressed on
    the waveform: complex gain of one tissue of a bulk time series (`sub`: one row per animal, gene counts and depth)."""
    g = list(wf.genes)
    cpm = sub[g].to_numpy(float) / sub.depth.to_numpy(float)[:, None] * 1e6
    ex = np.median(cpm, 0) >= 5
    c = cs.timed.cosinor(np.log(cpm[:, ex] + 0.5), zt)
    c = c["amplitude"] * np.exp(1j * W * c["acrophase_h"])
    x = (wf.amp * np.exp(1j * wf.psi))[ex]
    return complex(np.conj(x) @ c / (np.conj(x) @ x)), int(ex.sum())
