#!/usr/bin/env python3
"""
landauer_estimate.py
--------------------
The Landauer floor of the KSN variable and the timescale for approaching it
(Sections 2.5 and 2.6 of the paper, the KSN row of Table 1, and Section 3.2,
criterion 3).

Inputs: the energy and hashrate series in data/, and the quoted values below
(each with its source in the paper). The number of bits erased per hash
attempt, N_b, comes from the gate-level count in nb_gate_count.py.

Usage:  python3 landauer_estimate.py
"""

import numpy as np

from ksn_common import load_energy, annual_hashrate, HUBBLE_TIME_YR, T0

K_B = 1.380649e-23            # J/K
T_ROOM = 298.15               # K
GATES_PER_ATTEMPT = 240_896   # two-input gate operations per double SHA-256 attempt
N_B = 199_668                 # bits erased per attempt (nb_gate_count.py); "about 2e5"
BTC_TWH_PER_YR = 138.0        # network electricity use, 2024 [Neumueller2025]
FRONTIER = {                  # J per hash of the best available hardware
    2013: 600.0 / 66e9,       # first ASICs (Avalon: 66 GH/s at 600 W) [Taylor2017]
    2017: 0.07e-9,            # 16 nm chips, 0.07 W per GH/s [Taylor2017]
    2020: 29.5e-12,           # Antminer S19 Pro, 29.5 J/TH (range of rates only)
    2024: 13.5e-12,           # Antminer S21 XP, 13.5 J/TH [Bitmain2024]
}
KOOMEY_DOUBLING_YR = 2.7      # post-2000 doubling time [KoomeyNaffziger2015]
SIGNAL_KT = (30.0, 60.0)      # signal energy per operation for reliable switching [Agarwal2016]


def main():
    kT = K_B * T_ROOM
    print(f"k_B T = {kT:.3e} J;  k_B T ln 2 = {kT * np.log(2):.3e} J at {T_ROOM} K"
          "   [paper: 2.85e-21 J]")

    B_min = N_B * kT * np.log(2)
    print(f"N_b = {N_B} bits erased per hash attempt ({GATES_PER_ATTEMPT} gate operations)"
          "   [paper: N_b ~ 2e5]")
    print(f"B_min = N_b k_B T ln 2 = {B_min:.2e} J per hash, {np.log10(N_B):.1f} orders"
          " above the single-bit bound   [paper: ~6e-16; some five orders]")

    years_e, P = load_energy()
    years_h, H, _ = annual_hashrate()
    P_h = np.array([P[years_e == y][0] for y in years_h])
    B = P_h / H
    P24, H24, B24 = P_h[-1], H[-1], B[-1]
    f = BTC_TWH_PER_YR * 1e12 / (365.25 * 24) / P24
    eps = B24 * f
    print(f"\n2024: B = P/H = {B24:.3e} J per hash;  network {BTC_TWH_PER_YR:.0f} TWh/yr ->"
          f" f = {f:.1e} (1/f = {1 / f:.0f});  eps_hw = B f = {eps:.2e} J per hash")
    print("   [paper: f some 8e-4, 1/f ~ 10^3, eps_hw ~ 2.5e-11 J per hash]")
    print(f"gap B/B_min = {B24 / B_min:.1e} ({np.log10(B24 / B_min):.2f} dex);"
          f" hardware factor eps_hw/B_min = {eps / B_min:.1e};"
          f" allocation factor 1/f = {1 / f:.1e}")
    print("   [paper: seven to eight orders; 10^4 to 10^5; ~10^3]")

    print("\nDirect extrapolation of the trend of B in Figure 4 (unphysical):")
    for y0 in (2013, 2019):
        m = years_h >= y0
        slope, icpt = np.polyfit(years_h[m], np.log10(B[m]), 1)
        print(f"   fit {y0}-2024: {slope:+.3f} dex/yr -> B = B_min in "
              f"{(np.log10(B_min) - icpt) / slope:.1f}")
    H_need = P24 / B_min
    print(f"   needs H = P/B_min = {H_need:.1e} H/s ({H_need / H24:.1e} times 2024),"
          f" drawing {H_need * eps:.1e} W = {H_need * eps / P24:.1e} times world power")
    print("   [paper: between the 2040s and the 2060s; 10^4 to 10^5 times world power]")

    print("\nThe hardware factor: frontier energy per hash")
    for y, e in FRONTIER.items():
        print(f"   {y}: {e:.2e} J per hash")
    r_race = np.log10(FRONTIER[2013] / FRONTIER[2017]) / 4
    r_node = np.log10(FRONTIER[2017] / FRONTIER[2024]) / 7
    r_recent = np.log10(FRONTIER[2020] / FRONTIER[2024]) / 4
    r_koomey = np.log10(2) / KOOMEY_DOUBLING_YR
    print(f"   rates: 2013-2017 {r_race:.2f}, 2017-2024 {r_node:.3f}, 2020-2024 "
          f"{r_recent:.3f} dex/yr; Koomey-Naffziger {r_koomey:.3f} dex/yr")
    print("   [paper: about 1e-8 (2013), 7e-11 (2017), 1.35e-11 J (2024); about one"
          " order of magnitude per decade]")

    floors = [N_B * e * kT for e in SIGNAL_KT]
    print(f"\nPractical floor, {SIGNAL_KT[0]:.0f}-{SIGNAL_KT[1]:.0f} k_B T per operation:"
          f" {floors[0]:.1e} to {floors[1]:.1e} J per hash, {floors[0] / B_min:.0f} to"
          f" {floors[1] / B_min:.0f} times B_min ({np.log10(floors[0] / B_min):.2f} to"
          f" {np.log10(floors[1] / B_min):.2f} dex)")
    print("   [paper: some tens of k_B T; a few times 1e-14 J per hash; nearly two"
          " orders above B_min]")
    p_err = 0.01 / GATES_PER_ATTEMPT
    print(f"   (for comparison: one wrong hash in a hundred needs only"
          f" {np.log(1 / p_err):.0f} k_B T per operation)")

    drops = [np.log10(eps / fl) for fl in floors]          # dex still to fall
    t_fast = [d / r_koomey for d in drops]
    t_slow = [d / r_recent for d in drops]
    central_floor = N_B * np.mean(SIGNAL_KT) * kT
    t_central = np.log10(eps / central_floor) / r_node
    print(f"\nTime to the floor from 2024: {min(t_fast):.0f}-{max(t_fast):.0f} yr at"
          f" {r_koomey:.3f} dex/yr, {min(t_slow):.0f}-{max(t_slow):.0f} yr at"
          f" {r_recent:.3f} dex/yr; central (45 k_B T, {r_node:.3f} dex/yr):"
          f" {t_central:.0f} yr, year {2024 + t_central:.0f}")
    print("   [paper: about three decades; the 2050s]")
    plateau = [B24 / 10 ** d for d in drops]
    print(f"At fixed allocation, B falls by {min(drops):.1f}-{max(drops):.1f} dex to"
          f" {min(plateau):.1e}-{max(plateau):.1e} J per hash")
    print("   [paper: two and a half to three orders, to a few times 1e-11 J per hash]")
    dt = 2024 + t_central - T0
    print(f"\nTable 1, KSN row: 2050s; Delta t / (1/H0) = {dt / HUBBLE_TIME_YR:.1e}"
          "   [paper: ~6e-9]")


if __name__ == "__main__":
    main()
