#!/usr/bin/env python3
"""
convention_refit.py
-------------------
The energy-accounting conventions (Sections 2.1 and 3.3 of the paper).

The OWID column primary_energy_consumption follows the input-equivalent
(substitution) convention of the Energy Institute: electricity generated from
nuclear, hydro, wind and solar sources is counted as the fossil fuel that
would have been needed to generate it. This script
  1. shows that convention in the data (consumption / electricity ratios),
  2. rebuilds the series on the direct basis (that electricity counted as
     generated; fossil fuels, biofuels and other renewables at their input),
     and on the physical-content basis (nuclear electricity counted as the
     heat released in the reactor, at an efficiency of 33 percent),
  3. refits the linear and free-rate exponential models on each basis.

Usage:  python3 convention_refit.py
"""

import numpy as np
from scipy.optimize import curve_fit

from ksn_common import (load_energy_rows, TWH_PER_YR_TO_W, L_SUN, T0,
                        HUBBLE_TIME_YR)

NONCOMB = ("nuclear", "hydro", "wind", "solar")
NUCLEAR_THERMAL_EFFICIENCY = 0.33


def fits(P, t):
    """Linear (a, b, sigma_b, R^2) and exponential (r, sigma_r, R^2) fits."""
    b, a = np.polyfit(t, P, 1)
    res = P - (a + b * t)
    sst = np.sum((P - P.mean()) ** 2)
    sb = np.sqrt(np.sum(res ** 2) / (len(P) - 2) / np.sum((t - t.mean()) ** 2))
    (a0, r), cov = curve_fit(lambda tt, a0, r: a0 * np.exp(r * tt), t, P,
                             p0=(P[0], 0.02))
    R2e = 1 - np.sum((P - a0 * np.exp(r * t)) ** 2) / sst
    return a, b, sb, 1 - np.sum(res ** 2) / sst, r, np.sqrt(cov[1, 1]), R2e


def main():
    rows = load_energy_rows()
    years = np.array(sorted(rows))
    t = years - T0

    print("1. The published series is the substitution convention")
    print("   consumption / electricity generated, by source:")
    for y in (1965, 1990, 2000, 2010, 2024):
        r = rows[y]
        ratios = [f"{s} {r[s + '_consumption'] / r[s + '_electricity']:.2f}"
                  for s in NONCOMB if r[s + "_electricity"] > 0]
        print(f"   {y}: " + ", ".join(ratios))
    print("   (1/0.36 = 2.78 until 2000; about 2.4-2.5 by 2024, an efficiency of about"
          " 40 percent)   [paper: 36 percent until 2000, about 40 percent by 2024;"
          " some 2.5 to 2.8 times]")
    comp = np.array([rows[y]["fossil_fuel_consumption"] + rows[y]["nuclear_consumption"]
                     + rows[y]["renewables_consumption"] for y in years])
    pub = np.array([rows[y]["primary_energy_consumption"] for y in years])
    print(f"   check: fossil + nuclear + renewables = {100 * (comp / pub).min():.1f} to "
          f"{100 * (comp / pub).max():.1f} percent of the published total")

    def basis(nuclear_factor, other_ren_at="input"):
        out = []
        for y in years:
            r = rows[y]
            e = (r["fossil_fuel_consumption"] + r["biofuel_consumption"]
                 + (r["other_renewable_consumption"] if other_ren_at == "input"
                    else r["other_renewable_electricity"])
                 + nuclear_factor * r["nuclear_electricity"]
                 + r["hydro_electricity"] + r["wind_electricity"]
                 + r["solar_electricity"])
            out.append(e)
        return np.array(out)

    series = {
        "published (substitution)": pub,
        "direct": basis(1.0),
        "physical content": basis(1 / NUCLEAR_THERMAL_EFFICIENCY),
        "direct, other renewables at electricity (variant)": basis(1.0, "electricity"),
    }
    print("\n2. The series on each basis, relative to the published one")
    for name, E in series.items():
        if name.startswith("published"):
            continue
        print(f"   {name:50s} 1965: {100 * (E[0] / pub[0] - 1):+.1f}%   "
              f"2024: {100 * (E[-1] / pub[-1] - 1):+.1f}%")
    print("   [paper: on the direct basis 4 percent lower in 1965 and 11 percent lower in 2024]")
    phys, direct = series["physical content"], series["direct"]
    between = np.all((phys >= direct) & (phys <= pub))
    print(f"   physical content lies between direct and published in every year: {between}")

    print("\n3. Refits")
    for name, E in series.items():
        P = E * TWH_PER_YR_TO_W
        a, b, sb, R2l, r, sr, R2e = fits(P, t)
        tstar = (L_SUN - a) / b
        print(f"   {name:50s} b = ({b / 1e11:.2f} +/- {sb / 1e11:.2f})e11 W/yr, "
              f"R2 {R2l:.3f}; r = {100 * r:.2f} +/- {100 * sr:.2f} %/yr, R2 {R2e:.3f}; "
              f"Type II after {tstar / HUBBLE_TIME_YR:.2e} Hubble times")
    print("   [paper, direct basis: b = 2.16e11 W/yr, ~1.3e5 Hubble times, r = 1.94 %/yr;"
          " physical content between the two]")


if __name__ == "__main__":
    main()
