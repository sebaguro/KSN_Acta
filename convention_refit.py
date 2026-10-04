#!/usr/bin/env python3
"""
convention_refit.py
-------------------
The energy-accounting conventions (Sections 2.1 and 3.3 of the paper).

Energy statistics differ in how they count the electricity generated from
nuclear, hydro, wind and solar sources:

  substitution  as the fossil fuel that would have been needed to generate
                it. This is the convention of the series the paper analyses:
                the OWID column primary_energy_consumption in the release of
                27 April 2026 (data/owid_world_energy_1965_2024.csv).
  direct        as generated.
  third         nuclear electricity as the heat released in the reactor,
                and hydro, wind and solar electricity as generated. This is
                the convention of the total energy supply, which OWID has
                reported since its release of 10 September 2026
                (data/owid_world_total_energy_supply_1965_2024.csv).

This script
  1. shows the substitution convention in the series analysed
     (consumption / electricity ratios);
  2. builds the series under the other two conventions: the direct one by
     recounting the release analysed (the published total with the
     electricity of the four sources counted as generated, and nothing else
     changed), and the third one both as OWID publishes it and, for
     comparison, by recounting the release analysed (nuclear electricity
     counted as reactor heat at an efficiency of 33 percent);
  3. refits the linear and free-rate exponential models on each series by
     least squares, and compares the one-percent model with the linear fit
     as energy_fits.py does for the series analysed;
  4. repeats the model comparison of Section 2.4 on each series, with the
     sampler, seed and settings of waic_indistinguishability.py.

Usage:  python3 convention_refit.py
"""

import csv
import os

import numpy as np
from scipy.optimize import curve_fit

import waic_indistinguishability as waic
from ksn_common import (DATA, load_energy_rows, TWH_PER_YR_TO_W, L_SUN, T0,
                        HUBBLE_TIME_YR, KARDASHEV_RATE)

OWID_TES_CSV = os.path.join(DATA, "owid_world_total_energy_supply_1965_2024.csv")
NONCOMB = ("nuclear", "hydro", "wind", "solar")
NUCLEAR_THERMAL_EFFICIENCY = 0.33
TWH_YR_TO_TW = 1e12 / (365.25 * 24) / 1e12   # as in waic_indistinguishability.py


def fits(P, t):
    """Linear (a, b, sigma_b, R^2) and exponential (a0, r, sigma_r, R^2) fits,
    and the residual sum of squares of the one-percent model over that of the
    linear fit."""
    b, a = np.polyfit(t, P, 1)
    res = P - (a + b * t)
    sst = np.sum((P - P.mean()) ** 2)
    sb = np.sqrt(np.sum(res ** 2) / (len(P) - 2) / np.sum((t - t.mean()) ** 2))
    (a0, r), cov = curve_fit(lambda tt, a0, r: a0 * np.exp(r * tt), t, P,
                             p0=(P[0], 0.02))
    R2e = 1 - np.sum((P - a0 * np.exp(r * t)) ** 2) / sst
    # the one-percent model, anchored at the start of the record (Section 2.2)
    P0 = P[0] / (1 + KARDASHEV_RATE)
    ratio = np.sum((P - P0 * (1 + KARDASHEV_RATE) ** t) ** 2) / np.sum(res ** 2)
    return (a, b, sb, 1 - np.sum(res ** 2) / sst, a0, r, np.sqrt(cov[1, 1]), R2e,
            ratio)


def model_comparison(years, E):
    """Delta WAIC (exponential minus linear) and its standard error for the
    series E (TWh per year), computed as the unweighted analysis of
    waic_indistinguishability.py computes it: the same sampler, seed,
    starting points, step sizes and numbers of steps (its run_analysis())."""
    P = np.asarray(E, float) * TWH_YR_TO_TW          # TW
    t = (years - 1964).astype(float)
    rng = np.random.default_rng(waic.SEED)
    ll_lin, ll_exp = waic.make_loglikes(t, P, None)
    S_lin = waic.metropolis(rng, ll_lin, [5.0, 0.25, np.log(0.4)], waic.N_STEPS,
                            [0.15, 0.004, 0.08], waic.BURN)
    S_exp = waic.metropolis(rng, ll_exp, [5.0, 0.02, np.log(0.5)], waic.N_STEPS,
                            [0.15, 0.0007, 0.08], waic.BURN)
    e_lin = waic.pointwise_elpd(S_lin, t, P, "linear", None)
    e_exp = waic.pointwise_elpd(S_exp, t, P, "exponential", None)
    d = e_lin - e_exp
    dWAIC = (-2 * e_exp.sum()) - (-2 * e_lin.sum())   # >0: linear scored better
    se = 2.0 * np.sqrt(len(P)) * d.std(ddof=1)
    return dWAIC, se


def load_total_energy_supply(path=OWID_TES_CSV):
    """{year: {column: float}} for the World rows of the release of
    10 September 2026, 1965-2024."""
    rows = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows[int(r["year"])] = {k: float(v) for k, v in r.items()
                                    if k not in ("country", "year")}
    return rows


def span(values):
    return f"{min(values):.2f} to {max(values):.2f}"


def main():
    rows = load_energy_rows()
    years = np.array(sorted(rows))
    t = years - T0
    tes_rows = load_total_energy_supply()
    assert sorted(tes_rows) == list(years)

    print("1. The series analysed follows the substitution convention")
    print("   (OWID release of 27 April 2026, column primary_energy_consumption)")
    print("   consumption / electricity generated, by source:")
    for y in (1965, 1990, 2000, 2010, 2024):
        r = rows[y]
        ratios = [f"{s} {r[s + '_consumption'] / r[s + '_electricity']:.2f}"
                  for s in NONCOMB if r[s + "_electricity"] > 0]
        print(f"   {y}: " + ", ".join(ratios))
    together = [sum(rows[y][s + "_consumption"] for s in NONCOMB)
                / sum(rows[y][s + "_electricity"] for s in NONCOMB) for y in years]
    print(f"   the four sources together: {min(together):.3f} times the electricity generated"
          f" in {years[int(np.argmin(together))]}, the lowest value, and {max(together):.3f} in"
          f" {years[int(np.argmax(together))]}, the highest   [paper: some 2.4 to 2.8 times]")
    before = [100 / x for y, x in zip(years, together) if y < 2000]
    lo, hi = f"{min(before):.1f}", f"{max(before):.1f}"
    early = (f"{lo} percent in every year from 1965 to 1999" if lo == hi
             else f"{lo} to {hi} percent in the years 1965 to 1999")
    print(f"   as a conversion efficiency: {early},"
          f" {100 / together[list(years).index(2000)]:.1f} percent in 2000 and"
          f" {100 / together[-1]:.1f} percent in 2024"
          "   [paper: 36 percent until 2000, about 41 percent by 2024]")
    comp = np.array([rows[y]["fossil_fuel_consumption"] + rows[y]["nuclear_consumption"]
                     + rows[y]["renewables_consumption"] for y in years])
    pub = np.array([rows[y]["primary_energy_consumption"] for y in years])
    print(f"   check: fossil + nuclear + renewables = {100 * (comp / pub).min():.1f} to "
          f"{100 * (comp / pub).max():.1f} percent of the published total")
    print("   (the components do not add up exactly to the published total, so each recount"
          " below keeps that total and changes only the electricity of the four sources)")

    def recount(nuclear_factor, other_ren_at="input"):
        """The published total of the release analysed with the electricity of
        the four sources recounted: hydro, wind and solar as generated, nuclear
        as nuclear_factor times the electricity generated (1 under the direct
        convention, 1/0.33 for reactor heat)."""
        out = []
        for y in years:
            r = rows[y]
            e = r["primary_energy_consumption"]
            e -= r["nuclear_consumption"] - nuclear_factor * r["nuclear_electricity"]
            for src in ("hydro", "wind", "solar"):
                e -= r[src + "_consumption"] - r[src + "_electricity"]
            if other_ren_at == "electricity":
                e -= r["other_renewable_consumption"] - r["other_renewable_electricity"]
            out.append(e)
        return np.array(out)

    direct = recount(1.0)
    third_recount = recount(1 / NUCLEAR_THERMAL_EFFICIENCY)
    tes = np.array([tes_rows[y]["total_energy_twh"] for y in years])

    def against(E):
        return (f"1965: {100 * (E[0] / pub[0] - 1):+.1f}%   "
                f"2024: {100 * (E[-1] / pub[-1] - 1):+.1f}%")

    print("\n2. The other two conventions, relative to the series analysed")
    print("   a. direct: the release analysed, recounted")
    print(f"      {against(direct)}"
          "   [paper: 4 percent lower in 1965 and 10 percent lower in 2024]")
    print("      variant, other renewables also counted at their electricity:"
          f" {against(recount(1.0, 'electricity'))}")
    print("   b. third: the total energy supply of the OWID release of 10 September 2026,"
          " as published (column total_energy_twh)")
    print("      it follows the third convention: energy counted / electricity generated,"
          " 1965-2024:")
    for s in NONCOMB:
        ratios = [tes_rows[y][s + "_energy_twh"] / tes_rows[y][s + "_electricity_twh"]
                  for y in years if tes_rows[y][s + "_electricity_twh"] > 0]
        print(f"        {s:8s} {span(ratios)}   ({ratios[-1]:.2f} in 2024)")
    print(f"      {against(tes)}"
          "   [paper: 4 percent lower in 1965 and 7 percent lower in 2024]")
    print("   c. third: the release analysed, recounted (nuclear electricity as reactor heat,"
          f" efficiency {NUCLEAR_THERMAL_EFFICIENCY:.2f})")
    print(f"      {against(third_recount)}")
    gap = 100 * (1 - tes / third_recount)
    print(f"      the published total energy supply (b) is {gap.min():.2f} to {gap.max():.2f}"
          " percent lower than this recount")

    series = [
        ("substitution (the series analysed)", pub,
         "[paper: b = 2.44e11, ~1.1e5 Hubble times, r = 2.01]", "[paper: 4 +/- 18]"),
        ("direct (release analysed, recounted)", direct,
         "[paper: b = 2.18e11, ~1.3e5 Hubble times, r = 1.94]", "[paper: -4 +/- 16]"),
        ("third (total energy supply, as published)", tes,
         "[paper: b = 2.30e11, ~1.2e5 Hubble times, r = 1.96]",
         "[paper: 29 +/- 16; 1.8 standard errors]"),
        ("third (release analysed, recounted)", third_recount,
         "", "[paper: 27 +/- 17]"),
    ]

    print("\n3. Least-squares fits on each series")
    for name, E, quote, _ in series:
        P = E * TWH_PER_YR_TO_W
        a, b, sb, R2l, a0, r, sr, R2e, ratio = fits(P, t)
        tstar = (L_SUN - a) / b
        t_exp = T0 + np.log(L_SUN / a0) / r
        print(f"   {name}")
        print(f"      linear: b = ({b / 1e11:.2f} +/- {sb / 1e11:.2f})e11 W/yr, R2 {R2l:.3f};"
              f" Type II after {tstar / HUBBLE_TIME_YR:.2e} Hubble times")
        print(f"      exponential: r = {100 * r:.2f} +/- {100 * sr:.2f} %/yr, R2 {R2e:.3f};"
              f" Type II {(t_exp - 2024) / 1e3:.2f} kyr after 2024")
        print(f"      one-percent model: RSS(one-percent, anchored at 1965) / RSS(linear) ="
              f" {ratio:.0f}")
        if quote:
            print(f"      {quote}")
    print("   [paper: under each of the three conventions the one-percent rate is excluded]")

    print("\n4. Model comparison on each series (as in waic_indistinguishability.py, unweighted)")
    print("   Delta WAIC = WAIC(exponential) - WAIC(linear); positive values favour the linear form")
    largest = 0.0
    for name, E, _, quote in series:
        dWAIC, se = model_comparison(years, E)
        largest = max(largest, abs(dWAIC) / se)
        print(f"   {name:42s} Delta WAIC = {dWAIC:+6.2f}, standard error {se:5.2f},"
              f" ratio {dWAIC / se:+.2f}   {quote}")
    print("   (the first line repeats, for the series analysed, the unweighted result of"
          " waic_indistinguishability.py)")
    print(f"   largest |Delta WAIC| / standard error = {largest:.2f}: "
          + ("below 2 on every series, so neither form is significantly preferred under any"
             " convention" if largest < 2 else "2 or more on at least one series")
          + "   [paper: no significant preference under each of the three conventions]")


if __name__ == "__main__":
    main()
