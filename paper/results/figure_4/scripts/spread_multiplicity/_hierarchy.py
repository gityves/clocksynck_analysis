"""local: not in clocksync (paper-only model, REFACTOR_PLAN §1.3). Verbatim copy of legacy/src/sccc/hierarchy.py.

Hierarchical pooling of contrast estimates across diseases, cell families and cohorts (amendment A20).

For estimates y_i (one per contrast x compartment) with standard errors s_i:

    y_i ~ N(theta_i, s_i^2),  theta_i = mu_{D(i),F(i)} + a_{h(i)} + b_{c(i)} + e_i,  e_i ~ N(0, tau_e^2)
    mu_{D,F} ~ N(m_F, tau_D^2),  m_F ~ N(0, S^2),  a_h ~ N(0, tau_a^2) (shared control sets only),  b_c ~ N(0, tau_b^2)

with half-normal priors on the taus (scale `scale`; S = 2 scale). e_i is integrated out in the sampler
(y_i ~ N(mu + a + b, s_i^2 + tau_e^2)); theta_i is drawn afterwards from its conditional normal posterior, which gives
the local false sign rate lfsr_i = min(P(theta_i > 0), P(theta_i < 0)). `lfsr_specific` is the same for the
cell-specific part theta_i - a_h - b_c (post hoc, A20 run note): a shift of a whole contrast or cohort is not a
cell-type effect. `robust`: Student-t (nu = 4) likelihood on the
marginal residual (sensitivity; theta_i draws then use the normal approximation).
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd


def _setup(chains: int):
    os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1")
    import numpyro
    numpyro.set_host_device_count(chains)
    return numpyro


def _codes(x):
    lab, idx = np.unique(np.asarray(x).astype(str), return_inverse=True)
    return lab, idx


def fit(df: pd.DataFrame, y: str, se: str, scale: float, *, robust: bool = False, chains: int = 4,
        warmup: int = 1000, samples: int = 2000, seed: int = 0, shared=("MSSM", "HBCC", "RADC")) -> dict:
    """df columns: disease_class, family, contrast, cohort, `y`, `se`. Rows with a missing y or se are dropped.
    Returns {'cells': class x family summaries, 'rows': per-row theta summaries and lfsr, 'tau': variance
    components, 'diag': r_hat / divergences}."""
    numpyro = _setup(chains)
    import jax.numpy as jnp
    import numpyro.distributions as dist
    from jax import random
    from numpyro.infer import MCMC, NUTS

    d = df[np.isfinite(df[y]) & np.isfinite(df[se]) & (df[se] > 0)].reset_index(drop=True)
    cell_lab, cell = _codes(d.disease_class + "|" + d.family)
    fam_lab, fam = _codes(d.family)
    cell_fam = np.array([np.flatnonzero(fam_lab == c.split("|")[1])[0] for c in cell_lab])
    con_lab, con = _codes(d.contrast)
    coh = d.cohort.where(d.cohort.isin(shared), "__none__").to_numpy()
    coh_lab, coh_i = _codes(coh)
    has_coh = np.array([c != "__none__" for c in coh_lab])[coh_i].astype(float)
    yv, sv = d[y].to_numpy(float), d[se].to_numpy(float)

    def model():
        t_d = numpyro.sample("tau_D", dist.HalfNormal(scale))
        t_a = numpyro.sample("tau_a", dist.HalfNormal(scale))
        t_b = numpyro.sample("tau_b", dist.HalfNormal(scale))
        t_e = numpyro.sample("tau_e", dist.HalfNormal(scale))
        m = numpyro.sample("m_F", dist.Normal(0, 2 * scale).expand([len(fam_lab)]))
        zc = numpyro.sample("z_cell", dist.Normal(0, 1).expand([len(cell_lab)]))
        mu = numpyro.deterministic("mu", m[cell_fam] + t_d * zc)
        za = numpyro.sample("z_a", dist.Normal(0, 1).expand([len(coh_lab)]))
        zb = numpyro.sample("z_b", dist.Normal(0, 1).expand([len(con_lab)]))
        loc = mu[cell] + has_coh * t_a * za[coh_i] + t_b * zb[con]
        sd = jnp.sqrt(sv ** 2 + t_e ** 2)
        lik = dist.StudentT(4.0, loc, sd) if robust else dist.Normal(loc, sd)
        numpyro.sample("y", lik, obs=yv)
        numpyro.deterministic("loc", loc)

    mc = MCMC(NUTS(model, target_accept_prob=0.95), num_warmup=warmup, num_samples=samples, num_chains=chains,
              progress_bar=False, chain_method="parallel")
    mc.run(random.PRNGKey(seed))
    p = mc.get_samples()
    extra = mc.get_extra_fields()
    loc, te = np.asarray(p["loc"]), np.asarray(p["tau_e"])[:, None]
    prec = 1 / sv ** 2 + 1 / np.maximum(te, 1e-9) ** 2
    mean = (yv / sv ** 2 + loc / np.maximum(te, 1e-9) ** 2) / prec
    rng = np.random.default_rng(seed)
    theta = mean + rng.standard_normal(mean.shape) / np.sqrt(prec)
    pos = (theta > 0).mean(0)
    mu = np.asarray(p["mu"])
    spec = theta - (loc - mu[:, cell])  # cell-specific part mu + e (cohort and contrast effects removed)
    pos_s = (spec > 0).mean(0)
    rows = d.assign(theta_median=np.median(theta, 0), theta_lo=np.quantile(theta, 0.025, 0),
                    theta_hi=np.quantile(theta, 0.975, 0), p_pos=pos, lfsr=np.minimum(pos, 1 - pos),
                    specific_median=np.median(spec, 0), lfsr_specific=np.minimum(pos_s, 1 - pos_s),
                    shift_median=np.median(loc - mu[:, cell], 0))
    cells = pd.DataFrame({"cell": cell_lab, "disease_class": [c.split("|")[0] for c in cell_lab],
                          "family": [c.split("|")[1] for c in cell_lab], "median": np.median(mu, 0),
                          "lo": np.quantile(mu, 0.025, 0), "hi": np.quantile(mu, 0.975, 0), "p_pos": (mu > 0).mean(0),
                          "n_rows": np.bincount(cell, minlength=len(cell_lab)),
                          "n_contrasts": [d.contrast[cell == j].nunique() for j in range(len(cell_lab))]})
    tau = {k: float(np.median(np.asarray(p[k]))) for k in ("tau_D", "tau_a", "tau_b", "tau_e")}
    try:
        import arviz as az
        rhat = float(np.nanmax(az.rhat(az.from_numpyro(mc), var_names=["mu", "tau_D", "tau_e"]).to_array()))
    except Exception:  # diagnostics only
        rhat = np.nan
    diag = dict(divergences=int(np.asarray(extra["diverging"]).sum()), rhat_max=rhat, n=len(d))
    return dict(cells=cells, rows=rows, tau=tau, diag=diag)
