"""Unit checks of a220_lib against statsmodels PHReg (unpenalised Cox), a toy Harrell C and the lasso KKT conditions."""
import sys
from pathlib import Path

import numpy as np
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a220_lib as A  # noqa: E402

rng = np.random.default_rng(1); n, p = 400, 6
X = rng.normal(size=(n, p)); lp = X @ np.array([0.7, -0.5, 0.3, 0, 0, 0])
t = np.round(rng.exponential(np.exp(-lp)) * 10, 0) + 1  # rounded: many ties
d = (rng.random(n) < 0.7).astype(float); o = np.argsort(t, kind="stable"); X, t, d = X[o], t[o], d[o]
f, l = A.ties(t)
b, C, ok = A.cox_newton(X, d, f, l, 25, 0.0)
r = sm.PHReg(t, X, status=d, ties="breslow").fit()
print("newton vs PHReg max|dbeta|", np.abs(b - r.params).max(), "max|dse|", np.abs(np.sqrt(np.diag(C)) - r.bse).max())
assert np.abs(b - r.params).max() < 1e-5 and np.abs(np.sqrt(np.diag(C)) - r.bse).max() < 1e-5
pu = A.univariate_p(X, d, f, l); r1 = sm.PHReg(t, X[:, :1], status=d, ties="breslow").fit()
print("univariate p gene0", pu[0], r1.pvalues[0]); assert abs(pu[0] - r1.pvalues[0]) < 1e-6
# lasso: lambda -> 0 reproduces the unpenalised fit; KKT at an intermediate lambda
Xs = (X - X.mean(0)) / X.std(0); bs, _, _ = A.cox_newton(Xs, d, f, l, 25, 0.0)
B = A.lasso_path(Xs, d, f, l, np.array([1e-7]), 1e-12, 200)[0]; print("lasso(lam~0) vs newton", np.abs(B - bs).max()); assert np.abs(B - bs).max() < 1e-3
lm = A.lambda_max(Xs, d, f, l); print("lambda_max", lm, "beta at lmax*1.001", np.abs(A.lasso_path(Xs, d, f, l, np.array([lm * 1.001]))[0]).max())
lam = 0.3 * lm; Bl = A.lasso_path(Xs, d, f, l, np.array([lm, lam]), 1e-12, 200)[1]
ll, g, w = A.cox_stats(Xs @ Bl, d, f, l); gr = Xs.T @ g / n
act = Bl != 0; print("active", act.sum(), "KKT active", np.abs(gr[act] - lam * np.sign(Bl[act])).max() if act.any() else None, "KKT inactive max|grad|/lam", np.abs(gr[~act]).max() / lam)
assert np.abs(gr[act] - lam * np.sign(Bl[act])).max() < 1e-4 and np.abs(gr[~act]).max() <= lam * (1 + 1e-4)
# Harrell C toy: t=[1,2,3,4], d=[1,1,0,1], s=[4,3,2,1] -> pairs (0,1),(0,2),(0,3),(1,2),(1,3) all concordant; (3) no later; C=1; tie case
tt, dd = np.array([1., 2, 3, 4]), np.array([1., 1, 0, 1])
assert A.harrell_c(tt, dd, np.array([4., 3, 2, 1])) == 1.0 and A.harrell_c(tt, dd, np.array([1., 2, 3, 4])) == 0.0 and A.harrell_c(tt, dd, np.zeros(4)) == 0.5
assert abs(A.harrell_c(tt, dd, np.array([4., 3, 3, 1])) - 0.9) < 1e-12  # one tied comparable pair (1,2) counts 0.5 of 5
print("C toy ok; C(true lp) on simulated:", A.harrell_c(t, d, X @ np.array([0.7, -0.5, 0.3, 0, 0, 0])))
# log-rank sanity vs statsmodels
g = X[:, 0] > 0; print("logrank p", A.logrank_p(t, d, g), "statsmodels", sm.duration.survdiff(t, d, g.astype(int))[1])
print("ALL OK")
