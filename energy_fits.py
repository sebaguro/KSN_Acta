#!/usr/bin/env python3
"""
energy_fits.py
--------------
The analysis of the energy record (Sections 2.1-2.3, 2.7 and 3.1-3.3 of the
paper, and Tables 1 and 2), from the Our World in Data (OWID) World series
for 1965-2024 kept in data/.

Prints each quantity next to the value the paper states, as [paper: ...].
The Markov chain Monte Carlo (MCMC) posterior and the Widely Applicable
Information Criterion (WAIC) comparison are in waic_indistinguishability.py;
the energy-accounting conventions of Section 3.3 in convention_refit.py.

Usage:  python3 energy_fits.py
"""

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

from ksn_common import (load_energy, fractional_sigma, L_SUN, TYPE_I_1964,
                        TYPE_II_1964, T0, KARDASHEV_RATE, HUBBLE_TIME_YR,
                        AGE_UNIVERSE_YR, HABITABILITY_CEILING_W)

LINE = "-" * 72


def ols(t, P, w=None):
    """Linear fit P = a + b t (weights w = 1/sigma^2 if given).
    Returns a, b, formal standard error of b, R^2 and the residuals."""
    w = np.ones_like(P) if w is None else w
    X = np.vstack([np.ones_like(t), t]).T
    A = X.T @ (w[:, None] * X)
    a, b = np.linalg.solve(A, X.T @ (w * P))
    res = P - (a + b * t)
    R2 = 1 - np.sum(w * res ** 2) / np.sum(w * (P - np.average(P, weights=w)) ** 2)
    cov = np.linalg.inv(A) * (np.sum(res ** 2) / (len(P) - 2) if np.all(w == 1) else 1.0)
    return a, b, np.sqrt(cov[1, 1]), R2, res


def acf(x, k):
    x = x - x.mean()
    return np.sum(x[k:] * x[:-k]) / np.sum(x * x)


def ljung_box(x, m):
    n = len(x)
    Q = n * (n + 2) * sum(acf(x, k) ** 2 / (n - k) for k in range(1, m + 1))
    return Q, stats.chi2.sf(Q, m)


def durbin_watson(e):
    return np.sum(np.diff(e) ** 2) / np.sum(e ** 2)


def runs_test(x):
    """Wald-Wolfowitz runs test on the signs of x about its mean."""
    s = x > x.mean()
    n1, n2 = s.sum(), (~s).sum()
    R = 1 + np.sum(s[1:] != s[:-1])
    mu = 2 * n1 * n2 / (n1 + n2) + 1
    var = 2 * n1 * n2 * (2 * n1 * n2 - n1 - n2) / ((n1 + n2) ** 2 * (n1 + n2 - 1))
    z = (R - mu) / np.sqrt(var)
    return z, 2 * stats.norm.sf(abs(z))


def main():
    years, P = load_energy()
    t = years - T0
    N = len(P)

    print(LINE)
    print("Section 2.1  Data")
    print(LINE)
    print(f"N = {N} annual points, {years[0]}-{years[-1]}            [paper: N=60]")
    print(f"P(1965) = {P[0] / 1e12:.2f} TW                        [paper: 4.95 TW]")
    print(f"P(2024) = {P[-1] / 1e12:.2f} TW (= {P[-1] / 1e12:.1f} TW)       [paper: 20.2 TW]")

    print("\n" + LINE)
    print("Section 2.2  The one-percent model")
    print(LINE)
    g = (1 + KARDASHEV_RATE) ** t
    a, b, sb, R2l, res_l = ols(t, P)
    rss_lin = np.sum(res_l ** 2)
    P0 = P[0] / (1 + KARDASHEV_RATE)        # anchored at the start of the record
    rss_k = np.sum((P - P0 * g) ** 2)
    print(f"RSS(one-percent, anchored at 1965) / RSS(linear) = {rss_k / rss_lin:.0f}"
          "   [paper: more than two orders of magnitude]")
    P0f = np.sum(P * g) / np.sum(g * g)
    print(f"  (diagnostic: with its amplitude fitted instead, the ratio is "
          f"{np.sum((P - P0f * g) ** 2) / rss_lin:.0f})")
    t_k = np.log(TYPE_II_1964 / TYPE_I_1964) / KARDASHEV_RATE
    print(f"Kardashev's own inputs (4e12 W -> 4e26 W at 1%/yr, continuous): "
          f"{t_k:.0f} yr after 1964 -> year {T0 + t_k:.0f}")
    print("                                          [paper: ~3,200 yr; year ~5200]")
    R_E, S0 = 6.371e6, 1361.0
    print(f"Solar insolation at Earth, S0 pi R^2 = {S0 * np.pi * R_E ** 2:.3e} W"
          "        [paper: 1.74e17 W]")

    print("\n" + LINE)
    print("Section 2.3  The linear model")
    print(LINE)
    print(f"b = ({b / 1e11:.2f} +/- {sb / 1e11:.2f}) x 10^11 W/yr"
          f"  ({100 * sb / b:.2f} percent)   [paper: (2.44 +/- 0.04)e11, 1.5 percent]")
    print(f"a = {a:.3e} W,  R^2 = {R2l:.3f}                   [paper: R^2 = 0.987]")
    dP = np.diff(P)
    W, p = stats.shapiro(dP)
    print(f"Shapiro-Wilk on Delta P: W = {W:.3f}, p = {p:.4f}, skewness = "
          f"{stats.skew(dP):.3f}   [paper: 0.925, 0.0014, -0.664]")
    neg = years[1:][dP < 0]
    print(f"Years of negative growth: {', '.join(map(str, neg))}"
          "      [paper: e.g. 2009 and 2020]")
    for y in (2008, 2009, 2020):
        print(f"  Delta P({y}) = {dP[years[1:] == y][0] / 1e9:+.0f} GW")
    t_star = (L_SUN - a) / b
    print(f"Linear Type II time t* = (L_sun - a)/b = {t_star:.3e} yr"
          "        [paper: ~1.6e15 yr]")
    print(f"  in Hubble times (1/H0 = {HUBBLE_TIME_YR / 1e9:.2f} Gyr): "
          f"{t_star / HUBBLE_TIME_YR:.3e}      [paper: ~1.1e5]")
    print(f"  in ages of the Universe (13.8 Gyr): {t_star / AGE_UNIVERSE_YR:.3e}"
          "      [paper: five orders of magnitude]")

    print("\n" + LINE)
    print("Section 2.4  The free-rate exponential (least squares)")
    print(LINE)
    (a0, r), _ = curve_fit(lambda tt, a0, r: a0 * np.exp(r * tt), t, P,
                           p0=(P[0], 0.02))
    R2e = 1 - np.sum((P - a0 * np.exp(r * t)) ** 2) / np.sum((P - P.mean()) ** 2)
    print(f"r = {100 * r:.2f} %/yr, R^2 = {R2e:.3f}, Delta R^2 = {R2l - R2e:.4f}"
          "   [paper: 2.01, 0.987, < 0.001]")
    print("(the posterior r = 2.01 +/- 0.03 %/yr is in waic_indistinguishability.py)")
    r_post = 0.0201
    e = np.exp(r_post * t)
    a0p = np.sum(P * e) / np.sum(e * e)
    t2e = np.log(L_SUN / a0p) / r_post
    print(f"Type II year at r = 2.01 %/yr: {T0 + t2e:.0f}, i.e. "
          f"{(T0 + t2e - 2024) / 1e3:.2f} kyr after 2024   [paper: ~3500; some 1.5 kyr]")

    print("\n" + LINE)
    print("Section 2.7  Serial dependence of the growth rates xi = ln(P_t/P_t-1)")
    print(LINE)
    xi = np.diff(np.log(P))
    print(f"lag-one autocorrelation rho_1 = {acf(xi, 1):.2f}              [paper: 0.35]")
    for m in range(1, 5):
        Q, pq = ljung_box(xi, m)
        print(f"  Ljung-Box, lags 1-{m}: Q = {Q:.1f}, p = {pq:.4f}")
    print("                                   [paper: Q = 7.8-12.8, p < 0.02]")
    print(f"Durbin-Watson = {durbin_watson(xi - xi.mean()):.2f}"
          "                        [paper: 1.25]")
    z, pz = runs_test(xi)
    print(f"Runs test (about the mean): z = {z:.1f}, p = {pz:.3f}"
          "      [paper: z = -2.5, p = 0.014]")
    print("Diagnostics, not quoted in the paper: the same statistics once the trend"
          " in the growth rate is removed")
    tr = stats.linregress(years[1:], xi)
    print(f"  trend of xi: {100 * tr.slope:+.3f} percentage points per yr "
          f"(t = {tr.slope / tr.stderr:.1f}, p = {tr.pvalue:.4f}); fitted rate "
          f"{100 * (tr.intercept + tr.slope * 1966):.1f}% (1966) -> "
          f"{100 * (tr.intercept + tr.slope * 2024):.1f}% (2024)")
    e_xi = xi - (tr.intercept + tr.slope * years[1:])
    lb = [ljung_box(e_xi, m)[1] for m in range(1, 5)]
    zt, pzt = runs_test(e_xi)
    print(f"  after removing that trend: rho_1 = {acf(e_xi, 1):.2f}; Ljung-Box p "
          f"(lags 1-4) = {', '.join(f'{v:.3f}' for v in lb)}; Durbin-Watson = "
          f"{durbin_watson(e_xi):.2f}; runs p = {pzt:.3f}")
    lbd = [ljung_box(dP, m)[1] for m in range(1, 5)]
    print(f"  absolute increments Delta P: rho_1 = {acf(dP, 1):.2f}; Ljung-Box p "
          f"(lags 1-4) = {', '.join(f'{v:.2f}' for v in lbd)}; Durbin-Watson = "
          f"{durbin_watson(dP - dP.mean()):.2f}")

    print("\n" + LINE)
    print("Sections 2.4 and 3.1  The habitability ceiling (Balbi and Lingam 2025)")
    print(LINE)
    for c in HABITABILITY_CEILING_W:
        print(f"ceiling {c:.0e} W: {np.log10(L_SUN / c):.1f} orders below L_sun, "
              f"{np.log10(c / P[-1]):.1f} orders above P(2024)")
    print("                  [paper: eleven to twelve orders below L_sun;"
          " one to two orders of present power]")

    print("\n" + LINE)
    print("Section 3.3  Weighted re-analysis (adopted per-year uncertainties)")
    print(LINE)
    sig = fractional_sigma(years) * P
    w = 1 / sig ** 2
    aw, bw, _, R2w, res_w = ols(t, P, w)
    print(f"weighted b = {bw:.3e} W/yr ({100 * (bw / b - 1):+.1f} percent), "
          f"R^2 = {R2w:.3f}   [paper: 2.56e11, about 5 percent, 0.977]")
    chi_k = np.sum(w * (P - P0 * g) ** 2)
    print(f"weighted RSS(one-percent, anchored) / RSS(linear) = "
          f"{chi_k / np.sum(w * res_w ** 2):.0f}   [paper: more than two orders]")

    print("\n" + LINE)
    print("Tables 1 and 2  (timescales measured from 1964; 1/H0 = "
          f"{HUBBLE_TIME_YR / 1e9:.1f} Gyr)")
    print(LINE)
    rows = [("Kardashev 1964 (fixed r)", T0 + t_k, t_k),
            ("Exp. free r (r = 2.01%)", T0 + t2e, t2e),
            ("Linear OLS", T0 + t_star, t_star)]
    for name, year, dt in rows:
        print(f"{name:26s} Type II year {year:.3g}   Delta t / (1/H0) = "
              f"{dt / HUBBLE_TIME_YR:.1e}")
    print("   [paper, Table 1: ~5200, ~2e-7; ~3500, ~1e-7; ~1.6e15, ~1.1e5]")
    print("   (the KSN row of Table 1 is computed in landauer_estimate.py)")
    print(f"Table 2: Kardashev 1964 estimate 4 TW, Type II after {t_k / 1e3:.1f} kyr;"
          f" OWID 1965 {P[0] / 1e12:.2f} TW, {t_star / HUBBLE_TIME_YR:.1e} Hubble times")


if __name__ == "__main__":
    main()
