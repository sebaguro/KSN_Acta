#!/usr/bin/env python3
"""
rms_residuals.py

Root-mean-square (rms) residual of the linear and the free-rate exponential
fits to world primary energy, 1965-2024 (Sections 2.3 and 2.4 of the paper),
and the chi-squared and reduced chi-squared of the same fits against the
adopted per-year uncertainties of Section 3.3.

Put this file in the top folder of the KSN_Acta repository, next to
energy_fits.py, and run it from there:

    python3 rms_residuals.py

It reads one file, data/owid_world_energy_1965_2024.csv (the same data file
energy_fits.py uses; make_owid_extract.py rebuilds it from Our World in
Data's own file), writes nothing and uses no network. Needs numpy and scipy.
The fits are made the same way as in energy_fits.py. The paper quotes the
root-mean-square residuals in Sections 2.3 and 2.4 and in the Conclusions;
the chi-squared values are not quoted in the paper.
"""

import csv
import os

import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(HERE, "data", "owid_world_energy_1965_2024.csv")

TW = 1e12                                  # 1 TW in W
TWH_PER_YR_TO_W = 1e12 / (365.25 * 24.0)   # Eq. (1): TWh per year -> W
T0 = 1964                                  # time origin: t = year - 1964

# 1. Data: world primary energy E (TWh per year) converted to power P (W).
years, E = [], []
with open(CSV_FILE, newline="") as f:
    for row in csv.DictReader(f):
        years.append(int(row["year"]))
        E.append(float(row["primary_energy_consumption"]))
years = np.array(years)
P = np.array(E) * TWH_PER_YR_TO_W
t = years - T0
N = len(P)

# 2. Linear model by ordinary least squares: P_L(t) = a + b t.
b, a = np.polyfit(t, P, 1)
fit_lin = a + b * t

# 3. Free-rate exponential by least squares: P_E(t) = a0 exp(r t).
(a0, r), _ = curve_fit(lambda tt, a0, r: a0 * np.exp(r * tt), t, P,
                       p0=(P[0], 0.02))
fit_exp = a0 * np.exp(r * t)

# 4. Residuals, e_i = P_i - fit_i, in TW.
e_lin = (P - fit_lin) / TW
e_exp = (P - fit_exp) / TW

print("Year      P    linear  resid   resid^2    expon.  resid   resid^2")
print("          TW     TW     TW      TW^2        TW     TW      TW^2")
for i in range(N):
    print(f"{years[i]}  {P[i] / TW:6.3f}  {fit_lin[i] / TW:6.3f} "
          f"{e_lin[i]:+6.3f}  {e_lin[i] ** 2:7.4f}    "
          f"{fit_exp[i] / TW:6.3f} {e_exp[i]:+6.3f}  {e_exp[i] ** 2:7.4f}")

# 5. Root-mean-square residual: rms = sqrt( (1/N) * sum of e_i^2 ).
P_mean = P.mean() / TW
print()
print(f"N = {N} years ({years[0]}-{years[-1]}); mean power = {P_mean:.3f} TW")
print(f"Linear fit:      a = {a / TW:.4f} TW, b = {b:.4e} W/yr")
print(f"Exponential fit: a0 = {a0 / TW:.4f} TW, r = {100 * r:.3f} %/yr")
for name, e, paper in (("linear", e_lin, "0.48 TW, 4%"),
                       ("exponential", e_exp, "0.49 TW, 4%")):
    S = np.sum(e ** 2)
    rms = np.sqrt(S / N)
    print()
    print(f"{name}:")
    print(f"  sum of squared residuals   = {S:.4f} TW^2")
    print(f"  divided by N = {N}          = {S / N:.5f} TW^2")
    print(f"  square root = rms residual = {rms:.4f} TW")
    print(f"  rms / mean power           = {100 * rms / P_mean:.2f} %"
          f"      [paper: {paper}]")
    print(f"  with N - 2 = {N - 2} instead of N: s = {np.sqrt(S / (N - 2)):.4f} TW"
          "  (standard error of the regression)")

# 6. How R^2 relates to the rms residual. R^2 = 1 - sum(e^2) / sum((P - mean)^2),
#    so rms = sd(P) * sqrt(1 - R^2), where sd(P) is the spread of P about its
#    mean over the record. R^2 therefore gives the scatter only as a fraction
#    of that spread, not its size.
sdP = np.sqrt(np.mean((P / TW - P_mean) ** 2))
print()
print(f"spread of P about its mean, sd(P) = {sdP:.4f} TW")
for name, e in (("linear", e_lin), ("exponential", e_exp)):
    R2 = 1 - np.sum(e ** 2) / np.sum((P / TW - P_mean) ** 2)
    print(f"{name}: R^2 = {R2:.4f};  sd(P) * sqrt(1 - R^2) = "
          f"{sdP * np.sqrt(1 - R2):.4f} TW (= rms residual)")

# 7. Chi-squared against the adopted uncertainties of Section 3.3:
#    sigma_i = f_i P_i, with f_i falling linearly from 0.10 in 1965 to 0.015
#    in 2000 and constant thereafter.
#    chi^2 = sum (e_i / sigma_i)^2,  nu = N - p (p = 2),  chi^2_nu = chi^2 / nu.
f_i = np.clip(np.interp(years, [1965, 2000], [0.10, 0.015]), 0.015, 0.10)
sigma = f_i * P / TW
nu = N - 2
print()
print(f"Adopted uncertainties: sigma_i = f_i P_i; nu = N - 2 = {nu}")
periods = ((1965, 1979), (1980, 1999), (2000, 2024))
for name, e in (("linear", e_lin), ("exponential", e_exp)):
    z2 = (e / sigma) ** 2
    parts = ", ".join(f"{lo}-{hi}: {np.sum(z2[(years >= lo) & (years <= hi)]):.1f}"
                      for lo, hi in periods)
    print(f"{name}, unweighted fit: chi^2 = {np.sum(z2):.1f}, "
          f"chi^2_nu = {np.sum(z2) / nu:.2f}   (by period: {parts})")

# 8. The same with the weighted fits, which minimise chi^2 (weights 1/sigma_i^2).
b_w, a_w = np.polyfit(t, P / TW, 1, w=1 / sigma)
(a0_w, r_w), _ = curve_fit(lambda tt, a0, r: a0 * np.exp(r * tt), t, P / TW,
                           p0=(P[0] / TW, 0.02), sigma=sigma, absolute_sigma=True)
chi2_lin_w = np.sum(((P / TW - (a_w + b_w * t)) / sigma) ** 2)
chi2_exp_w = np.sum(((P / TW - a0_w * np.exp(r_w * t)) / sigma) ** 2)
print(f"linear, weighted fit (b = {b_w * TW:.3e} W/yr): chi^2 = {chi2_lin_w:.1f}, "
      f"chi^2_nu = {chi2_lin_w / nu:.2f}, sqrt(chi^2_nu) = {np.sqrt(chi2_lin_w / nu):.2f}")
print(f"exponential, weighted fit (r = {100 * r_w:.2f} %/yr): chi^2 = {chi2_exp_w:.1f}, "
      f"chi^2_nu = {chi2_exp_w / nu:.2f}, sqrt(chi^2_nu) = {np.sqrt(chi2_exp_w / nu):.2f}")
print(f"weighted fits: chi^2_exp - chi^2_lin = {chi2_exp_w - chi2_lin_w:.1f}"
      "   (compare Delta WAIC = -50 +/- 38, Section 3.3)")
print(f"expected scatter of chi^2_nu for a correct model: sqrt(2/nu) = {np.sqrt(2 / nu):.2f}")
