#!/usr/bin/env python3
"""
intermediate_forms.py
---------------------
Growth laws between the linear and the exponential (response to Reviewer 2;
not fitted in the paper, decision D7). Fits a power law and a stretched
exponential to the energy record alongside the linear and free-rate
exponential forms, and gives for each the year it would reach the waste-heat
habitability ceiling of Balbi and Lingam (2025), of order 3e14 to 1e15 W, and
the Type II threshold, L_sun.

It supports the statement in Section 2.4: a trajectory growing faster than
linearly but slower than exponentially crosses both thresholds at
intermediate times and in the same order.

Usage:  python3 intermediate_forms.py
"""

import numpy as np
from scipy.optimize import brentq, curve_fit

from ksn_common import load_energy, L_SUN, T0, HABITABILITY_CEILING_W


def linear(t, a, b):
    return a + b * t


def exponential(t, a0, r):
    return a0 * np.exp(r * t)


def power_law(t, A, t0, alpha):
    return A * (1 + t / t0) ** alpha


def stretched(t, A, tau, beta):
    return A * np.exp((t / tau) ** beta)


MODELS = [
    ("linear", linear, lambda P: (P[0], 2.4e11)),
    ("power law", power_law, lambda P: (P[0], 60.0, 2.0)),
    ("stretched exponential", stretched, lambda P: (P[0] / 2, 30.0, 0.7)),
    ("exponential", exponential, lambda P: (P[0], 0.02)),
]
NAMES = {"power law": "A, t0, alpha", "stretched exponential": "A, tau, beta",
         "linear": "a, b", "exponential": "a0, r"}


def crossing(f, p, level, t_start=60.0):
    """Time t (years after 1964) at which f(t) first reaches level, t > t_start."""
    hi = t_start + 1.0
    while f(hi, *p) < level:
        hi *= 2
        if hi > 1e20:
            return np.inf
    return brentq(lambda t: f(t, *p) - level, t_start, hi, xtol=1e-6, rtol=1e-12)


def main():
    years, P = load_energy()
    t = years - T0
    sst = np.sum((P - P.mean()) ** 2)
    print(f"{'form':22s} {'R^2':>7s}  parameters")
    fits = {}
    for name, f, p0 in MODELS:
        p, _ = curve_fit(f, t, P, p0=p0(P), maxfev=200000)
        fits[name] = (f, p)
        R2 = 1 - np.sum((P - f(t, *p)) ** 2) / sst
        print(f"{name:22s} {R2:7.4f}  {NAMES[name]} = " + ", ".join(f"{v:.3g}" for v in p))

    print("\nYears after 2024 at which each form reaches")
    head = "".join(f"{c:>14.0e} W" for c in HABITABILITY_CEILING_W)
    print(f"{'form':22s}{head}   Type II (L_sun)")
    for name, (f, p) in fits.items():
        c = [crossing(f, p, lev) + T0 - 2024 for lev in HABITABILITY_CEILING_W]
        t2 = crossing(f, p, L_SUN) + T0 - 2024
        print(f"{name:22s}" + "".join(f"{x:16.3g}" for x in c) + f"   {t2:.3g}")
    print("\nEvery form reaches the habitability ceiling long before the Type II"
          " threshold; the growth law changes the times, not the order.")


if __name__ == "__main__":
    main()
