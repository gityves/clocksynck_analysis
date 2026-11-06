"""A220 survival kernels (numba): Breslow Cox partial likelihood, unpenalised Newton fit with Wald p, L1-penalised Cox by cyclic
coordinate descent (glmnet-style, standardised columns, objective -loglik/n + lambda |beta|_1), Harrell's C, two-group log-rank.
All Cox routines assume observations SORTED by time ascending; `ties(t)` gives the first / last index of each tie group."""
import math

import numba as nb
import numpy as np


@nb.njit(cache=True)
def ties(t):
    n = len(t); first = np.empty(n, np.int64); last = np.empty(n, np.int64); i = 0
    while i < n:
        j = i
        while j + 1 < n and t[j + 1] == t[i]:
            j += 1
        for k in range(i, j + 1):
            first[k] = i; last[k] = j
        i = j + 1
    return first, last


@nb.njit(cache=True)
def cox_stats(eta, d, first, last):
    """Breslow log partial likelihood, gradient and diagonal Hessian weights w (>= 0) with respect to eta."""
    n = len(eta); m = eta.max(); r = np.exp(eta - m)
    S = np.empty(n); acc = 0.0
    for i in range(n - 1, -1, -1):
        acc += r[i]; S[i] = acc
    ca = np.zeros(n); cb = np.zeros(n); a = 0.0; b = 0.0; ll = 0.0
    for i in range(n):
        s = S[first[i]]
        if d[i] > 0:
            ll += eta[i] - m - np.log(s); a += 1.0 / s; b += 1.0 / (s * s)
        ca[i] = a; cb[i] = b
    g = np.empty(n); w = np.empty(n)
    for i in range(n):
        A = ca[last[i]]; B = cb[last[i]]
        g[i] = d[i] - r[i] * A; w[i] = r[i] * A - r[i] * r[i] * B
    return ll, g, w


@nb.njit(cache=True)
def cox_ll(eta, d, first, last):
    return cox_stats(eta, d, first, last)[0]


@nb.njit(cache=True)
def cox_newton(X, d, first, last, iters=25, ridge=1e-6):
    """Unpenalised Cox (tiny ridge for stability): (beta, covariance of beta, converged)."""
    n, p = X.shape; beta = np.zeros(p); ok = False
    cur = cox_ll(X @ beta, d, first, last)
    H = np.eye(p)
    for _ in range(iters):
        eta = X @ beta; m = eta.max(); r = np.exp(eta - m)
        grad = np.zeros(p); H = np.zeros((p, p)); S0 = 0.0; S1 = np.zeros(p); S2 = np.zeros((p, p)); i = n - 1
        while i >= 0:
            f = first[i]; dg = 0.0; xs = np.zeros(p)
            for k in range(f, i + 1):
                S0 += r[k]; S1 += r[k] * X[k]; S2 += r[k] * np.outer(X[k], X[k])
                if d[k] > 0:
                    dg += 1.0; xs += X[k]
            if dg > 0:
                grad += xs - dg * S1 / S0; H += dg * (S2 / S0 - np.outer(S1, S1) / (S0 * S0))
            i = f - 1
        Hr = H + ridge * np.eye(p); step = np.linalg.solve(Hr, grad - ridge * beta)
        t = 1.0; moved = False
        for _ in range(12):
            nb_ = beta + t * step; v = cox_ll(X @ nb_, d, first, last) - 0.5 * ridge * np.sum(nb_ * nb_)
            if v >= cur - 1e-12:
                beta = nb_; cur = v; moved = True; break
            t *= 0.5
        if not moved or np.max(np.abs(t * step)) < 1e-8:
            ok = moved or np.max(np.abs(step)) < 1e-6; break
    return beta, np.linalg.inv(H + ridge * np.eye(p)), ok


@nb.njit(cache=True)
def univariate_p(X, d, first, last):
    """Wald p of each column in a one-covariate Cox model."""
    p = X.shape[1]; out = np.ones(p)
    for j in range(p):
        b, C, ok = cox_newton(X[:, j:j + 1].copy(), d, first, last, 25, 1e-8)
        if ok and C[0, 0] > 0:
            out[j] = math.erfc(abs(b[0]) / math.sqrt(C[0, 0]) / math.sqrt(2.0))
    return out


@nb.njit(cache=True)
def lasso_path(X, d, first, last, lams, tol=1e-7, maxouter=40):
    """L1 Cox path (warm starts) on standardised X (columns mean 0, variance 1): minimises -loglik/n + lam |beta|_1. Rows = lams."""
    n, p = X.shape; Bs = np.zeros((len(lams), p)); beta = np.zeros(p); eta = np.zeros(n)
    for l in range(len(lams)):
        lam = lams[l]
        for _ in range(maxouter):
            ll, g, w = cox_stats(eta, d, first, last)
            for i in range(n):
                w[i] = max(w[i], 1e-9)
            den = np.zeros(p)
            for j in range(p):
                s = 0.0
                for i in range(n):
                    s += w[i] * X[i, j] * X[i, j]
                den[j] = s / n
            res = g / w; old = beta.copy(); oldobj = -ll / n + lam * np.abs(beta).sum()
            for _ in range(200):
                dmax = 0.0
                for j in range(p):
                    num = 0.0
                    for i in range(n):
                        num += w[i] * X[i, j] * (res[i] + X[i, j] * beta[j])
                    num /= n
                    z = (abs(num) - lam) if abs(num) > lam else 0.0
                    bn = math.copysign(z, num) / den[j]; db = bn - beta[j]
                    if db != 0.0:
                        for i in range(n):
                            res[i] -= X[i, j] * db
                        beta[j] = bn
                        if den[j] * db * db > dmax:
                            dmax = den[j] * db * db
                if dmax < tol * 0.1:
                    break
            t = 1.0; dl = beta - old
            for _ in range(8):  # step halving if the penalised objective went up
                nb_ = old + t * dl; e2 = X @ nb_
                if -cox_ll(e2, d, first, last) / n + lam * np.abs(nb_).sum() <= oldobj + 1e-12:
                    break
                t *= 0.5
            beta = old + t * dl; eta = X @ beta
            if np.max(den * (beta - old) ** 2) < tol:
                break
        Bs[l] = beta
    return Bs


@nb.njit(cache=True)
def lambda_max(X, d, first, last):
    n = X.shape[0]; ll, g, w = cox_stats(np.zeros(n), d, first, last)
    return np.abs(X.T @ g).max() / n


@nb.njit(cache=True)
def harrell_c(t, d, s):
    """Harrell's C: comparable pairs (i event, t_j > t_i); concordant if s_i > s_j, half for tied scores. NaN with no comparable pair."""
    conc = 0.0; tot = 0.0
    for i in range(len(t)):
        if d[i] > 0:
            for j in range(len(t)):
                if t[j] > t[i]:
                    tot += 1.0
                    if s[i] > s[j]:
                        conc += 1.0
                    elif s[i] == s[j]:
                        conc += 0.5
    return conc / tot if tot > 0 else np.nan


def logrank_p(t, d, g):
    """Two-group log-rank chi-square p (g boolean)."""
    from scipy.stats import chi2
    o = e = v = 0.0
    for tt in np.unique(t[d > 0]):
        at = t >= tt; n = at.sum(); n1 = (at & g).sum(); dd = ((t == tt) & (d > 0)).sum(); d1 = ((t == tt) & (d > 0) & g).sum()
        o += d1; e += dd * n1 / n
        if n > 1:
            v += dd * (n1 / n) * (1 - n1 / n) * (n - dd) / (n - 1)
    return float(chi2.sf((o - e) ** 2 / v, 1)) if v > 0 else 1.0


def standardise(Xtr, Xte=None):
    mu, sd = Xtr.mean(0), Xtr.std(0); keep = sd > 1e-12; mu, sd = mu[keep], sd[keep]
    out = [(Xtr[:, keep] - mu) / sd]
    if Xte is not None:
        out.append((Xte[:, keep] - mu) / sd)
    return keep, out


def lasso_cv(X, t, d, folds, K=25, ratio=0.01):
    """L1 Cox on training data (sorted by t; X standardised): lambda by CV partial likelihood (glmnet's cross-validated deviance:
    held-out contribution ll(all; beta_-k) - ll(rest; beta_-k), summed over folds, maximised). Returns (beta at lambda.min, lambda index, lams)."""
    f, l = ties(t)
    lmax = lambda_max(X, d, f, l)
    if not lmax > 0:
        return np.zeros(X.shape[1]), 0, np.array([0.0])
    lams = lmax * np.exp(np.linspace(0, np.log(ratio), K))
    full = lasso_path(X, d, f, l, lams)
    score = np.zeros(K)
    for k in np.unique(folds):
        tr = folds != k; ft, lt = ties(t[tr]); Bk = lasso_path(X[tr], d[tr], ft, lt, lams)
        for j in range(K):
            eta = X @ Bk[j]
            score[j] += cox_ll(eta, d, f, l) - cox_ll(eta[tr], d[tr], ft, lt)
    j = int(np.argmax(score))
    return full[j], j, lams
