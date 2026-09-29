#!/usr/bin/env python3
"""
waic_indistinguishability.py
----------------------------
Demonstrates that a linear OLS model and a free-rate exponential model are
statistically indistinguishable descriptions of the 1965-2024 global
primary-energy record, by computing WAIC for both models and -- the crucial
step -- the standard error of the WAIC difference.

Two modes, matching the paper:

  (default)    UNWEIGHTED  -- the primary analysis of Section 2.4. A single
               homoscedastic noise scale sigma is sampled as a nuisance
               parameter (equivalent to plain OLS for the trend parameters).
               This is the assumption-free baseline: OWID publishes no
               per-year uncertainties.

  --weighted   WEIGHTED    -- the Section 3.3 sensitivity check. Adopts the
               paper's heteroscedastic error model (fractional 1-sigma
               declining from 10% in 1965 to 1.5% from 2000 on,
               sigma_i = f_i * P_i, as in make_errorbar_figures.py) and
               refits both models under fixed 1/sigma_i^2 weighting.

  --both       run both modes one after the other.

Expected behaviour (the point of the exercise):
  * Unweighted: Delta WAIC = about +4 (linear side), SE ~ 18  ->  ~0.2 SE.
  * Weighted:   Delta WAIC ~ -50 (exponential side), SE ~ 38  ->  ~1.3 SE.
  The mild, non-significant edge FLIPS SIGN with the error treatment --
  fractional-error weighting is nearly a log-space fit, the exponential's
  home turf -- while neither treatment ever approaches the 2*SE threshold.
  The only conclusion robust to the error model is indistinguishability.

Method references:
  WAIC              : Watanabe (2010), J. Mach. Learn. Res. 11, 3571
  SE of Delta-WAIC  : Vehtari, Gelman & Gabry (2017), Stat. Comput. 27, 1413

Data: OWID energy dataset (Ritchie et al.), column `primary_energy_consumption`
      (input-equivalent, or substitution, convention of the Energy Institute),
      country == "World", 1965-2024 (N = 60); the copy in data/ (see
      data/README.md for its version).

Sampler: Metropolis-Hastings, 75,000 steps of which the first 15,000 are
      discarded as burn-in, keeping 60,000 posterior samples, as stated in
      Section 2.4 of the paper; prior r ~ N(0.01, 0.02^2); random seed 7.

Usage:
  python3 waic_indistinguishability.py [--weighted | --both] [csv]
"""

import argparse
import csv
import os
import numpy as np
from scipy.special import logsumexp

SEED = 7
N_STEPS = 75_000          # 60,000 samples kept after the burn-in
BURN = 15_000
YEAR_RANGE = (1965, 2024)

# ---------------------------------------------------------------------------
# 1. Data loading: TWh/yr -> TW  (Eq. 1 of the paper)
# ---------------------------------------------------------------------------
def load_world_energy(path):
    twh_yr_to_tw = 1e12 / (365.25 * 24) / 1e12   # TWh/yr -> W -> TW
    rows = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r.get("country", "").strip() == "World":
                try:
                    y = int(float(r["year"]))
                    v = r["primary_energy_consumption"]
                    if v not in ("", None):
                        rows[y] = float(v) * twh_yr_to_tw          # TW
                except (KeyError, ValueError):
                    pass
    yrs = np.array(sorted(y for y in rows if YEAR_RANGE[0] <= y <= YEAR_RANGE[1]))
    P = np.array([rows[y] for y in yrs])
    t = (yrs - 1964).astype(float)       # paper's time origin: t = year - 1964
    return yrs, t, P

# ---------------------------------------------------------------------------
# 2. The paper's adopted heteroscedastic error model (Section 3.3;
#    identical to fractional_sigma() in make_errorbar_figures.py).
#    These are ADOPTED sensitivity-test uncertainties, not measurements.
# ---------------------------------------------------------------------------
def fractional_sigma(years):
    hi, lo, y0, y1 = 0.10, 0.015, 1965, 2000
    return np.clip(np.interp(years, [y0, y1], [hi, lo]), lo, hi)

# ---------------------------------------------------------------------------
# 3. Closed-form weighted linear OLS (validates against Sec 3.3 numbers)
# ---------------------------------------------------------------------------
def weighted_ols(t, P, sig):
    w = 1.0 / sig**2
    X = np.vstack([np.ones_like(t), t]).T
    W = np.diag(w)
    a, b = np.linalg.inv(X.T @ W @ X) @ (X.T @ W @ P)
    res = P - (a + b * t)
    R2w = 1.0 - np.sum(w * res**2) / np.sum(w * (P - np.average(P, weights=w))**2)
    return a, b, R2w

# ---------------------------------------------------------------------------
# 4. Metropolis-Hastings MCMC
# ---------------------------------------------------------------------------
def metropolis(rng, loglike, p0, steps, scales, burn, thin=1):
    p = np.array(p0, float)
    lp = loglike(p)
    keep = []
    for i in range(steps):
        q = p + np.array(scales) * rng.standard_normal(len(p))
        lq = loglike(q)
        if np.log(rng.random()) < lq - lp:
            p, lp = q, lq
        if i >= burn and i % thin == 0:
            keep.append(p.copy())
    return np.array(keep)

# ---------------------------------------------------------------------------
# 5. Log-likelihoods.
#    Unweighted: parameters (trend..., log sigma); sigma sampled, same for
#                every year -> equivalent to OLS for the trend.
#    Weighted:   parameters (trend...) only; sigma_i FIXED to the adopted
#                model, identically for BOTH models (symmetric treatment).
# ---------------------------------------------------------------------------
def make_loglikes(t, P, sig=None):
    if sig is None:                                          # UNWEIGHTED
        def ll_linear(theta):
            a, b, ls = theta
            s = np.exp(ls)
            mu = a + b * t
            return np.sum(-0.5 * np.log(2*np.pi*s*s) - 0.5*((P-mu)/s)**2)

        def ll_expon(theta):
            a0, r, ls = theta
            if a0 <= 0:
                return -1e18
            s = np.exp(ls)
            mu = a0 * np.exp(r * t)
            prior_r = -0.5 * ((r - 0.01) / 0.02) ** 2        # r ~ N(0.01, 0.02^2)
            return np.sum(-0.5*np.log(2*np.pi*s*s) - 0.5*((P-mu)/s)**2) + prior_r
    else:                                                    # WEIGHTED
        const = -0.5 * np.log(2 * np.pi * sig * sig)

        def ll_linear(theta):
            a, b = theta
            return np.sum(const - 0.5 * ((P - (a + b*t)) / sig) ** 2)

        def ll_expon(theta):
            a0, r = theta
            if a0 <= 0:
                return -1e18
            prior_r = -0.5 * ((r - 0.01) / 0.02) ** 2
            return np.sum(const - 0.5*((P - a0*np.exp(r*t)) / sig)**2) + prior_r
    return ll_linear, ll_expon

# ---------------------------------------------------------------------------
# 6. WAIC from pointwise predictive densities:
#      lppd_i   = log-mean-exp_s  log p(y_i | theta_s)
#      p_waic_i = Var_s  log p(y_i | theta_s)
#      elpd_i   = lppd_i - p_waic_i ;   WAIC = -2 * sum_i elpd_i
# ---------------------------------------------------------------------------
def pointwise_elpd(samples, t, P, model, sig=None):
    S, N = len(samples), len(P)
    L = np.empty((S, N))
    for j, th in enumerate(samples):
        if sig is None:
            *trend, ls = th
            s = np.exp(ls)
        else:
            trend, s = th, sig
        if model == "linear":
            a, b = trend
            mu = a + b * t
        else:
            a0, r = trend
            mu = a0 * np.exp(r * t)
        L[j] = -0.5 * np.log(2*np.pi*s*s) - 0.5 * ((P - mu) / s) ** 2
    lppd = logsumexp(L, axis=0) - np.log(S)
    p_waic = L.var(axis=0)
    return lppd - p_waic

# ---------------------------------------------------------------------------
# 7. One full analysis (either mode)
# ---------------------------------------------------------------------------
def run_analysis(yrs, t, P, weighted):
    rng = np.random.default_rng(SEED)            # re-seed per mode: each mode
    N = len(P)                                   # is independently reproducible
    label = "WEIGHTED (Sec 3.3 sensitivity)" if weighted else "UNWEIGHTED (Sec 2.4 primary)"
    print("=" * 66)
    print(f"  {label}")
    print("=" * 66)

    sig = fractional_sigma(yrs) * P if weighted else None

    if weighted:
        a_w, b_w, R2w = weighted_ols(t, P, sig)
        print(f"Weighted linear OLS: b = {b_w*1e12:.3e} W/yr  [paper Sec 3.3: 2.56e11]")
        print(f"                     weighted R^2 = {R2w:.3f}   [paper Sec 3.3: 0.977]\n")

    ll_lin, ll_exp = make_loglikes(t, P, sig)
    if weighted:
        S_lin = metropolis(rng, ll_lin, [5.0, 0.25], N_STEPS, [0.05, 0.003], BURN)
        S_exp = metropolis(rng, ll_exp, [5.0, 0.02], N_STEPS, [0.05, 0.0006], BURN)
    else:
        S_lin = metropolis(rng, ll_lin, [5.0, 0.25, np.log(0.4)], N_STEPS,
                           [0.15, 0.004, 0.08], BURN)
        S_exp = metropolis(rng, ll_exp, [5.0, 0.02, np.log(0.5)], N_STEPS,
                           [0.15, 0.0007, 0.08], BURN)

    r_m, r_s = S_exp[:, 1].mean(), S_exp[:, 1].std()
    b_m = S_lin[:, 1].mean()
    anchor = "[paper: 2.01 +/- 0.03]" if not weighted else ""
    print(f"Posterior r (exponential) = {100*r_m:.2f}% +/- {100*r_s:.2f}% per yr  {anchor}")
    lo, hi = np.percentile(S_exp[:, 1], [2.5, 97.5])
    print(f"  95% credible interval   = [{100*lo:.2f}%, {100*hi:.2f}%]"
          f"  {'[paper: [1.94%, 2.08%]]' if not weighted else ''}")
    print(f"  samples kept            = {len(S_exp)}")
    print(f"Posterior b (linear)      = {b_m*1e12:.3e} W/yr"
          f"  {'[paper: 2.44e11]' if not weighted else ''}\n")

    e_lin = pointwise_elpd(S_lin, t, P, "linear", sig)
    e_exp = pointwise_elpd(S_exp, t, P, "exponential", sig)
    WAIC_lin, WAIC_exp = -2 * e_lin.sum(), -2 * e_exp.sum()
    d = e_lin - e_exp
    dWAIC = WAIC_exp - WAIC_lin                  # >0: linear scored better
    se = 2.0 * np.sqrt(N) * d.std(ddof=1)        # Vehtari et al. (2017), x2 for deviance

    print(f"WAIC_linear      = {WAIC_lin:9.2f}")
    print(f"WAIC_exponential = {WAIC_exp:9.2f}")
    print(f"Delta WAIC (exp - lin) = {dWAIC:+8.2f}"
          f"   {'[paper: about 4]' if not weighted else '[paper Sec 3.3: -50]'}")
    print(f"SE(Delta WAIC)         = {se:8.2f}"
          f"   {'[paper: ~18]' if not weighted else '[paper Sec 3.3: 38]'}")
    print(f"Delta WAIC / SE        = {dWAIC/se:+8.2f}")
    verdict = ("STATISTICALLY INDISTINGUISHABLE (|Delta WAIC| < 2 SE)"
               if abs(dWAIC) < 2 * se else "genuine model preference (>2 SE)")
    side = "linear" if dWAIC > 0 else "exponential"
    print(f"\n=> {verdict}; the non-significant edge sits on the {side} side.\n")

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="?",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "data", "owid_world_energy_1965_2024.csv"),
                    help="path to the OWID energy CSV (default: the copy in data/)")
    ap.add_argument("--weighted", action="store_true",
                    help="run the Sec 3.3 weighted sensitivity analysis")
    ap.add_argument("--both", action="store_true",
                    help="run unweighted then weighted")
    args = ap.parse_args()

    yrs, t, P = load_world_energy(args.csv)
    print(f"Loaded N = {len(P)} points, {yrs[0]}-{yrs[-1]}, "
          f"P(1965) = {P[0]:.2f} TW, P(2024) = {P[-1]:.2f} TW\n")

    if args.both:
        run_analysis(yrs, t, P, weighted=False)
        run_analysis(yrs, t, P, weighted=True)
    else:
        run_analysis(yrs, t, P, weighted=args.weighted)
